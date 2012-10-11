#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "highOrder.h"
#include "lowOrder.h"
#include "fct.h"
#include "reg_estimator.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = ((M_c)_ij/tau + theta D_ij) (u_i^H - u_j^H) - ((M_c)_ij/tau - (1-theta) D_ij) (u_i^n -u_j^n) ; f^new_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  (M_L/tau -theta(K+D)) u^(n+1) = M_L/tau u^L + f^new 


const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX = 5; 
                     
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.
 
const double EPS_smooth = 1e-5;   				//for smoothness-indicator :  a < b =>  a+EPS < b
const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Adaptivity
const double THRESHOLD_UNREF = 0.001; 							// Unrefinement: error of all sons is smaller than THRESHOLD_UNREF times maximum element error
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_METHOD = 1;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.3;                      // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hanging nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 2.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 5;												// max. numbers of adaptivity steps

//Visualization
const bool HERMES_VISUALIZATION = false;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 1000;													//Every VTK_FREQth time step the solution is saved as VTK output.


// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";



int main(int argc, char* argv[])
{  

   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);


  // Initialize boundary conditions.
  DefaultEssentialBCConst<double>  bc_essential(BDY_IN, 0.0);
  EssentialBCs<double>  bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
	H1Space<double> space(&mesh, &bcs, P_INIT);	

 // Initialize solutions  
  Solution<double>  low_sln, ref_sln, high_sln, sln;
  

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step);
	CustomWeakFormConvection  convection;

  // Initialize views.
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	OrderView mview("mesh", new WinGeom(500, 0, 500, 400));
	OrderView ref_mview("ref_mesh", new WinGeom(0, 0, 500, 400));
	char title[100];
	
  // Output solution in VTK format.
	Linearizer lin;
	Orderizer ord;
	bool mode_3D = true;


  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

	//Initialize
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	double* u_L = NULL; 
	double* u_H =NULL;
	double* ref_sln_double =NULL;
	

	int ref_ndof, ndof; double err_est_rel_total;
	Adapt<double>* adaptivity = new Adapt<double>(&space, HERMES_L2_NORM);
	OGProjection<double> ogProjection;	
	Lumped_Projection lumpedProjection;	
	Low_Order lowOrder(theta);
	High_Order highOrd(theta);
	Flux_Correction fluxCorrection(theta);
	Regularity_Estimator regEst(EPS_smooth);
	
	DiscreteProblem<double>* dp_mass = new DiscreteProblem<double>(&massmatrix, &space);
	DiscreteProblem<double>* dp_convection = new DiscreteProblem<double>(&convection, &space); 

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	
do
{
	 info("Current Time: %f, Time step %d",current_time, ts);	 
// Spatial adaptivity loop. Note: u_prev_time must not be changed during spatial adaptivity. 
    // Periodic global derefinement. 
   if ((ts > 1 && ts % UNREF_FREQ == 0)||(space.get_num_dofs() >= NDOF_STOP)) 
    { 
    	info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: Exceptions::Exception("Wrong global derefinement method.");
      }      
    }
 
    bool done = false; int as = 1;	
	
			
    do 
    {
      	info("Time step %d, adaptivity step %d:", ts, as);	
				
				//Unrefinement
				if(as==1){
					ndof = space.get_num_dofs(); 
					double* coeff_vec_smooth = new double[ndof];
					OGProjection<double>::project_global(&space,&u_prev_time, coeff_vec_smooth, matrix_solver, HERMES_L2_NORM);				  
					Solution<double>::vector_to_solution(coeff_vec_smooth, &space, &low_sln);
   				err_est_rel_total = adaptivity->calc_err_est(&low_sln, &u_prev_time) * 100;    
  				adaptivity->unrefine(THRESHOLD_UNREF);		
  				delete [] coeff_vec_smooth; 
				  // Visualize the solution.
  				if(HERMES_VISUALIZATION)
  				{
						sprintf(title, "unrefined Mesh: Time %3.2f,timestep %i", current_time,ts);
						mview.set_title(title); 
						mview.show(&space);
					}
    		}	
				ndof = space.get_num_dofs(); 
								   	    
			ndof = space.get_num_dofs();  	    		
			double* coeff_vec_smooth = new double[ndof];
			int* smooth_elem_ref;	
							
			//smoothness-check for projected data		
			ogProjection.project_global(&space,&u_prev_time, coeff_vec_smooth,matrix_solver, HERMES_L2_NORM);					
			smooth_elem_ref =regEst.get_smooth_elems(&space,coeff_vec_smooth);
			// Construct reference mesh and setup reference space.					
			Mesh*  ref_mesh = new Mesh;
			ref_mesh->copy(space.get_mesh());
			Space<double>* ref_space = space.dup(ref_mesh, 0);
			HPAdapt * adapting = new HPAdapt(ref_space, HERMES_L2_NORM);	
							// increase p in smooth regions, h refine in non-smooth regions 
			if(adapting->adapt_smooth(smooth_elem_ref, P_MAX)==false) 
				throw Exceptions::Exception("reference space couldn't be constructed");							
									      
			delete adapting;
			delete [] coeff_vec_smooth; 							
   

				ref_ndof = ref_space->get_num_dofs();
				
			if(HERMES_VISUALIZATION) 
			{
				sprintf(title, "Ref_Mesh: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
				ref_mview.set_title(title);
				ref_mview.show(ref_space);
			}
	
	double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];
			double* limited_flux = new double[ref_ndof];	
	
			dp_mass->set_spaces(ref_space);
			dp_convection->set_spaces(ref_space);
			fluxCorrection.init(ref_space);				

			dp_mass->assemble(mass_matrix); 										//M_c/tau
			dp_convection->assemble(conv_matrix, NULL,true);		//K

		//----------------------MassLumping  & Artificial Diffusion --------------------------------------------------------------------	
			UMFPackMatrix<double>* lumped_matrix = fluxCorrection.massLumping(mass_matrix); // M_L/tau
			UMFPackMatrix<double>* diffusion = fluxCorrection.artificialDiffusion(conv_matrix);	
						
			//-----------------Assembling of matrices ---------------------------------------------------------------------	
			lowOrder.assemble_Low_Order(conv_matrix,diffusion,lumped_matrix);	
			highOrd.assemble_High_Order(conv_matrix,mass_matrix);
	
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L


			//--------- Project the previous timestep solution on the FE space (FCT is applied )----------------			
			// coeff_vec : FCT -Projection, coeff_vec_2: L2 Projection (ogProjection)			
				fluxCorrection.project_FCT(&u_prev_time, coeff_vec, coeff_vec_2,mass_matrix,lumped_matrix,time_step,&ogProjection,&lumpedProjection, &regEst);			
	//------------------------- lower order solution------------					
			u_L = lowOrder.solve_Low_Order(lumped_matrix, coeff_vec,time_step);								
	//-------------high order solution (standard galerkin) ------				
			u_H = highOrd.solve_High_Order(coeff_vec);		
		//------------------------------Assemble antidiffusive fluxes and limit these-----------------------------------	
			fluxCorrection.antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, u_L,coeff_vec, limited_flux,time_step,&regEst);
	//-------------Compute final solution ---------------			
				ref_sln_double = lowOrder.explicit_Correction(limited_flux);
			Solution<double> ::vector_to_solution(ref_sln_double, ref_space, &ref_sln);	
			
			


	      // Project the fine mesh solution onto the coarse mesh.
	     OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_L2_NORM); 

	      // Calculate element errors and total error estimate.
	    err_est_rel_total = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

	      // Report results.
	      info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
		   ndof, ref_ndof, err_est_rel_total);

	      // If err_est too large, adapt the mesh.
	   	if (err_est_rel_total < ERR_STOP){
				//info("Adaptivity-Stop: err_est_rel_total");
			 	done = true;
			
	    }else if(as>=ADAPSTEP_MAX){
					//info("Adaptivity-Stop: ADAPSTEP_MAX");
			 		done = true;		
	    }else{ 	      
				//info("Adapting the coarse mesh.");
				done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);	
				//info("Adapted.");	
				if (Space<double>::get_num_dofs(&space) >= NDOF_STOP){ 
					info("Adaptivity-Stop:NDOF_STOP,ndof: %d,refdof: %d ",Space<double>::get_num_dofs(&space),Space<double>::get_num_dofs(ref_space));
					done = true;
				}else
					// Increase the counter of performed adaptivity steps.
					as++;
	    }
		    if(done==true)    // Copy last reference solution into u_prev_time. 
     u_prev_time.copy(&ref_sln);
	   
	      // Visualize the solution and mesh.
	  if(HERMES_VISUALIZATION)
	  {
				sprintf(title, "Ref-Loesung: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
				sview.set_title(title);
				sview.show(&ref_sln);
				sprintf(title, "Mesh: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
				mview.set_title(title);
				mview.show(&space);
			}
			
		 if((VTK_VISUALIZATION) &&((done==true)&&(ts  % VTK_FREQ == 0)))
		 {

			// Output solution in VTK format.
			  char filename[40];
			  sprintf(filename, "solution-%i.vtk", ts );
			  lin.save_solution_vtk(&u_prev_time, filename, "solution", mode_3D);  
			  sprintf(filename, "ref_space_order-%i.vtk", ts);
				ord.save_orders_vtk(ref_space, filename);
					sprintf(filename, "ref_mesh-%i.vtk", ts );
					ord.save_mesh_vtk(ref_space, filename);       
					
   	 } 
		

     
     
     if((done==true)&&(current_time +time_step >= T_FINAL)){
     CustomInitialCondition exact_solution(ref_space->get_mesh());
	Adapt<double>* error_estimation = new Adapt<double>(ref_space, HERMES_L2_NORM);	
double err_est = error_estimation->calc_err_est(&exact_solution,&u_prev_time,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
double err_est_2 = error_estimation->calc_err_est(&u_prev_time,&exact_solution,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
FILE * pFile;
pFile = fopen ("error.txt","w");
     fprintf (pFile, "err_est = %f, err_est_2 =%f,  ndof = %d", err_est,err_est_2, ref_ndof);
fclose (pFile);
     
     		lin.save_solution_vtk(&u_prev_time, "solution_end.vtk", "solution", mode_3D);
		ord.save_orders_vtk(ref_space, "ref_space_end.vtk");
		ord.save_mesh_vtk(ref_space, "ref_mesh_end.vtk");     
     delete error_estimation;
     
     }
		
	      // Clean up.
			delete lumped_matrix; 
			delete diffusion;
			delete [] coeff_vec_2;
			delete [] coeff_vec; 
			delete [] limited_flux; 

			delete ref_mesh; 
			delete ref_space; 
	  
	


	  
    }
    while (done == false);

  // Update global time.
  current_time += time_step;


  // Increase time step counter
  ts++;
}
while (current_time < T_FINAL); 



   // Save solution in VTK format.
		lin.save_solution_vtk(&u_prev_time, "solution_end.vtk", "solution", mode_3D);
		ord.save_orders_vtk(&space, "space_end.vtk");
		ord.save_mesh_vtk(&space, "mesh_end.vtk");

// Clean up.
	delete adaptivity;
	delete mass_matrix;  
	delete conv_matrix;
	delete dp_convection;
	delete dp_mass;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

