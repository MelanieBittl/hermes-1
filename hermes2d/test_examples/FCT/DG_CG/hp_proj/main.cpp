#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "fct.h"
#include "reg_estimator.h"
#include "prev_solution.h"
#include <list>



using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;



const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX = 2; 										//Maximal polynomial degree.           

const bool serendipity = true;
 
const double EPS_smooth = 1e-10;   		//constant for the smoothness indicator (a<b => a+eps<=b)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Adaptivity
 
const double THRESHOLD = 0.2;                      // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const CandList CAND_LIST = H2D_H_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 600000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 15;												// max. numbers of adaptivity steps


//Visualization
const bool HERMES_VISUALIZATION = false;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION =true;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 70000000;													//Every VTK_FREQth time step the solution is saved as VTK output.



#include "error_estimations.cpp"

int main(int argc, char* argv[])
{  
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", basemesh);
 //mloader.load("tri.mesh", basemesh);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);
  
  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT,serendipity));	

 // Initialize solution of lower & higher order
  MeshFunctionSharedPtr<double> ref_sln(new Solution<double>),sln(new Solution<double>);

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>exact_solution(new CustomInitialCondition(mesh, true));
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix;



  // Output solution in VTK format.
	Linearizer lin;
	Orderizer ord;
	bool mode_3D = true;
	char filename[40];
	char title[100];
//Hermes Visualization
				OrderView ref_mview("ref_mesh", new WinGeom(500, 0, 500, 400));
			  	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
			  	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));


	//Initialize
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;   //M_c/tau

	double* ref_sln_double =NULL;
	int ref_ndof, ndof; double err_est_rel_total;
	
	  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST,H2DRS_DEFAULT_ORDER);
       selector.set_error_weights(1.0,1.0,1.0); 
	
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
  
  Adapt<double> adaptivity(space, &error_calculator);
  //AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);
  AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
  adaptivity.set_strategy(&stoppingCriterion);

	OGProjection<double> ogProjection;	
	Lumped_Projection lumpedProjection;	
	Flux_Correction fluxCorrection;
	Regularity_Estimator regEst(EPS_smooth);

	DiscreteProblem<double> dp_mass(&massmatrix, space);



 
    bool done = false; int as = 1;

	
    do 
    	{			ndof = space->get_num_dofs();  
		  		Hermes::Mixins::Loggable::Static::info("adap_step %i, dof = %i,", as, ndof);				
							
 	
			    		
		double* coeff_vec_smooth = new double[ndof];
			int* smooth_elem_ref;	
							

//      Hermes::Mixins::Loggable::Static::info("Projecting...");
				ogProjection.project_global(space,exact_solution, coeff_vec_smooth, HERMES_L2_NORM);		

			
   Hermes::Mixins::Loggable::Static::info("Calling get_smooth_elems()...");		
			smooth_elem_ref =regEst.get_smooth_elems(space,coeff_vec_smooth);

      // Construct reference mesh and setup reference space->
	     		 MeshSharedPtr ref_mesh(new Mesh);
      ref_mesh->copy(space->get_mesh()); 
		//ref_mesh->refine_all_elements();  ///for h-adapt only
      Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh, 0);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();


   HPAdapt* adapting = new HPAdapt(ref_space);	
							// increase p in smooth regions, h refine in non-smooth regions 
    //  Hermes::Mixins::Loggable::Static::info("Calling adapt_smooth()...");
			if(adapting->adapt_smooth(smooth_elem_ref, P_MAX)==false) 
				throw Exceptions::Exception("reference space couldn't be constructed");							
									      
			delete adapting;
			delete [] coeff_vec_smooth; 	
			    

			ref_ndof = ref_space->get_num_dofs();
			
			if(HERMES_VISUALIZATION) 
			{
				sprintf(title, "Ref_Mesh: as=%i,",as);
				ref_mview.set_title(title);
				ref_mview.show(ref_space);
				//View::wait();
			}


 Hermes::Mixins::Loggable::Static::info("Ref space assigning:%i",ref_ndof);	
      dp_mass.set_space(ref_space);

			fluxCorrection.init(ref_space);	

			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];				

			dp_mass.assemble(mass_matrix); 										//M_c
		
		//----------------------MassLumping  & Artificial Diffusion --------------------------------------------------------------------	
			 // Hermes::Mixins::Loggable::Static::info("Mass lumping and art. diffusion assembling...");	
			CSCMatrix<double>* lumped_matrix = fluxCorrection.massLumping(mass_matrix); // M_L


			//--------- Project 	
			// coeff_vec : FCT -Projection, coeff_vec_2: L2 Projection (ogProjection)	

				fluxCorrection.project_FCT(exact_solution, coeff_vec, coeff_vec_2,mass_matrix,lumped_matrix,&ogProjection,&lumpedProjection, &regEst);			

			Solution<double>::vector_to_solution(coeff_vec, ref_space, ref_sln);			
			
			
		    // Project the fine mesh solution onto the coarse mesh->
		   ogProjection.project_global(space, ref_sln, sln, HERMES_L2_NORM); 
		    // Calculate element errors and total error estimate.
   	error_calculator.calculate_errors(sln, ref_sln);
		  err_est_rel_total = error_calculator.get_total_error_squared() * 100;
			    // Report results.
	 		Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", ndof,ref_ndof, err_est_rel_total);				
				  // If err_est_rel too large, adapt the mesh->
		  if((err_est_rel_total < ERR_STOP)||(as>=ADAPSTEP_MAX)) done = true;
		  else
		  {
		    done = adaptivity.adapt(&selector);
		    // Increase the counter of performed adaptivity steps.
		    if(done == false)  as++;
		  }
		  if(space->get_num_dofs() >= NDOF_STOP) done = true;	 
		 
					


    	 
	      // Visualize the solution and mesh->
	  if(HERMES_VISUALIZATION)
	  {

				sprintf(title, "Ref-Loesung: as=%i,", as);
				sview.set_title(title);
				sview.show(ref_sln);
				sprintf(title, "Mesh: as=%i,", as);
				mview.set_title(title);
				mview.show(space);
			}
			

   	 
   	if((VTK_VISUALIZATION)&&(done==true))
   	{
			calc_error_total(ref_sln, exact_solution,ref_space);
					ord.save_mesh_vtk(ref_space, "end_ref_mesh.vtk");
		ord.save_orders_vtk(ref_space, "end_ref_order.vtk");	
		lin.save_solution_vtk(ref_sln, "end_ref_solution2d.vtk", "solution", false);
   	lin.save_solution_vtk(ref_sln, "end_ref_solution3d.vtk", "solution", mode_3D);
	
   	
   	}
		
	      // Clean up.
			delete lumped_matrix; 
			delete [] coeff_vec_2;
			delete [] coeff_vec; 




	  
    }
    while (done == false);




	delete mass_matrix;  


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

