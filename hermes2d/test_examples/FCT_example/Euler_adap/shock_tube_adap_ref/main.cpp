#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "lumped_projection.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;



const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX =1;

const double time_step = 1e-4;
const double T_FINAL = 0.231;                       // Time interval length. 

const double theta = 0.5;

// Equation parameters.   
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;  

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


// Adaptivity
const double THRESHOLD_UNREF = 0.001; 							// Unrefinement: error of all sons is smaller than THRESHOLD_UNREF times maximum element error
const int UNREF_FREQ = 100;                         // Every UNREF_FREQth time step the mesh is derefined.
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
const CandList CAND_LIST = H2D_H_ANISO;          // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
const int ADAPSTEP_MAX = 5;												// max. numbers of adaptivity steps



//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"


//Visualization
const bool HERMES_VISUALIZATION = false;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool HERMES_VISUALIZATION_LOW = false;           // Set to "false" to suppress Hermes OpenGL visualization for low_order.
const bool HERMES_VISUALIZATION_ORDER = false;           // Set to "false" to suppress Hermes OpenGL visualization for low_order.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 500;													//Every VTK_FREQth time step the solution is saved as VTK output. 
     

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

int main(int argc, char* argv[])
{
   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain2.mesh", &basemesh);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

	 double diag;			
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8.;
	//info("diag=%f, h_min = %f", diag, h_min);

	H1Space<double> space_rho(&mesh,P_INIT);
  H1Space<double> space_rho_v_x(&mesh,P_INIT);
  H1Space<double> space_rho_v_y(&mesh,P_INIT);
  H1Space<double> space_e(&mesh, P_INIT);

	int dof_rho = space_rho.get_num_dofs();
	int dof_v_x = space_rho_v_x.get_num_dofs();
	int dof_v_y = space_rho_v_y.get_num_dofs();
	int dof_e = space_e.get_num_dofs();


  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
	CustomInitialCondition_rho prev_rho(&mesh);
  ConstantSolution<double> prev_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e prev_e(&mesh, KAPPA);

 Solution<double>  new_rho,new_rho_v_x,new_rho_v_y,new_e;
 Solution<double> coarse_rho,coarse_rho_v_x,coarse_rho_v_y,coarse_rho_e;
 
  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  
//InitialCondition auch fuer low setzen wegen Filtern!
	CustomInitialCondition_rho low_rho(&mesh);
  ConstantSolution<double> low_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> low_rho_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e low_rho_e(&mesh, KAPPA);

	CustomInitialCondition_rho boundary_rho(&mesh);
  ConstantSolution<double> boundary_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> boundary_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e boundary_e(&mesh,KAPPA);

		Solution<double> R_h_1, R_h_2;

 //--------- Filters for visualization of pressure & velocity
 	PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 2);

  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  s1.set_min_max_range(0, 1.);
  s2.set_min_max_range(0., 1.);
  s3.set_min_max_range(0., 1.);
	pressure_view.set_min_max_range(0.,1.);	
 /* if(HERMES_VISUALIZATION){  
			s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
			pressure_view.show(&pressure);
	}*/

 
	OrderView mview_r("mesh_rho", new WinGeom(0, 0, 500, 400));
	OrderView mview_vx("mesh_vx", new WinGeom(300, 0, 500, 400));
	OrderView mview_vy("mesh_vy", new WinGeom(0, 300, 500, 400));
	OrderView mview_e("mesh_e", new WinGeom(300, 300, 500, 400));


   Linearizer lin_p, lin_v_x, lin_v_y, lin_rho;
   Orderizer ord_p;		
  		
	PressureFilter pressure_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), KAPPA);
	VelocityFilter vel_x_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 1);
	VelocityFilter vel_y_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 2);
		ScalarView s1_n("low_rho", new WinGeom(0, 0, 600, 300));
		ScalarView s2_n("low_rho_v_x", new WinGeom(700, 0, 600, 300));
		ScalarView s3_n("low_rho_v_y", new WinGeom(0, 400, 600, 300));
		ScalarView s4_n("low_pressure", new WinGeom(700, 400, 600, 300));
     s1_n.set_min_max_range(0., 1.);
      s2_n.set_min_max_range(0., 1.);
      s3_n.set_min_max_range(0., 1.);
      s4_n.set_min_max_range(0., 1.);
//------------


//-----------Weakforms
  EulerEquationsWeakForm_Mass wf_mass(time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, time_step, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);
  EulerBoundary wf_boundary(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerBoundary wf_boundary_low(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);


  // Set up the solver, matrix, and rhs according to the solver selection.
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   
	UMFPackMatrix<double> * matrix_L_low = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_L = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS_low = new UMFPackMatrix<double> ; 

	double* u_L = NULL; 
	char title[100]; char filename[40];	
	AsmList<double>*  al = new AsmList<double>;	
	
// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	int as = 1;
	int ref_ndof =0;

 Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

do
{	 
		info(" Time step %d, time %3.5f, ndof=%i", ts, current_time, ndof); 
		
		
	if ((ts > 1 && ts % UNREF_FREQ == 0)||(space_rho.get_num_dofs() >= NDOF_STOP)) 
    { 
    	info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        case 2: mesh.unrefine_all_elements();
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        case 3: mesh.unrefine_all_elements();
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        default: Exceptions::Exception("Wrong global derefinement method.");
      }      
    }		
  
    bool done = false; as = 1;

	//Adaptivity loop
	do{
   	info("Time step %d, adaptivity step %d:", ts, as);	
   	
      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space<double> *>* ref_spaces 
          = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),0,0);
      Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1],(*ref_spaces)[2], (*ref_spaces)[3]);
   	
		//ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
		ref_ndof = Space<double>::get_num_dofs(ref_spaces_const);
		
		dof_rho = (*ref_spaces)[0]->get_num_dofs();
		dof_v_x = (*ref_spaces)[1]->get_num_dofs();
		dof_v_y = (*ref_spaces)[2]->get_num_dofs();
		dof_e = (*ref_spaces)[3]->get_num_dofs();

  // Initialize the FE problem.
  	DiscreteProblem<double>* dp_mass = new DiscreteProblem<double>(&wf_mass, ref_spaces_const);
 	 DiscreteProblem<double> *dp_boundary = new DiscreteProblem<double>(&wf_boundary, ref_spaces_const);
 	 DiscreteProblem<double>* dp_K = new DiscreteProblem<double>(&wf_K, ref_spaces_const);


			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ref_ndof);
			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];	
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];	


    dp_mass->assemble(mass_matrix);
	//----------------------MassLumping M_L--------------------------------------------------------------------
		UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of the initial condition
			Lumped_Projection::project_lumped(ref_spaces_const,Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec, matrix_solver);
			
	if(HERMES_VISUALIZATION_LOW){
		Solution<double>::vector_to_solutions(coeff_vec,ref_spaces_const, Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));	
			s1_n.show(&low_rho);
			s2_n.show(&vel_x_low);
			s3_n.show(&vel_y_low);
			s4_n.show(&pressure_low);
}					
			OGProjection<double>::project_global(ref_spaces_const,Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);


		dp_boundary->assemble(matrix_dS);
    dp_K->assemble(lowmat_rhs);

					//------------------------artificial DIFFUSION D---------------------------------------		
			UMFPackMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,(*ref_spaces)[0], (*ref_spaces)[1],(*ref_spaces)[2], (*ref_spaces)[3], dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);
			lowmat_rhs->add_matrix(diffusion); //L(U)
			lowmat_rhs->add_matrix(matrix_dS); //L(U)+dS(U) 

			matrix_L->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)
			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)

	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, ref_spaces_const, Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));
			 }else error("Matrix solver failed.\n");


		if(HERMES_VISUALIZATION_LOW){
				s1_n.show(&low_rho);
				s2_n.show(&vel_x_low);
				s3_n.show(&vel_y_low);
				s4_n.show(&pressure_low);
			}



		 	 DiscreteProblem<double>* dp_boundary_low = new DiscreteProblem<double>(&wf_boundary_low,ref_spaces_const);
				DiscreteProblem<double>* dp_K_low = new DiscreteProblem<double>(&wf_K_low,ref_spaces_const);

					dp_boundary_low->assemble(matrix_dS_low);	
			 	 	dp_K_low->assemble(matrix_L_low);
					UMFPackMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,(*ref_spaces)[0], (*ref_spaces)[1],(*ref_spaces)[2], (*ref_spaces)[3], dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
					matrix_L_low->add_matrix(diffusion_low); //L(U)
					matrix_L_low->add_matrix(matrix_dS_low); //L(U)+dS(U) 


				//---------------------------------------antidiffusive fluxes-----------------------------------	
					antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);
					for(int i=0; i<ref_ndof;i++){
									 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		
					}

					Solution<double>::vector_to_solutions(coeff_vec, ref_spaces_const, Hermes::vector<Solution<double> *>(&new_rho, &new_rho_v_x, &new_rho_v_y, &new_e));	


      // Project the fine mesh solution onto the coarse meshes.
info("Projecting fine mesh solutions on coarse meshes for error estimation.");
      	OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&new_rho, &new_rho_v_x, &new_rho_v_y, &new_e), 
	  Hermes::vector<Solution<double> *>(&coarse_rho,&coarse_rho_v_x,&coarse_rho_v_y,&coarse_rho_e)); 

   // Calculate element errors and total error estimate.
							info("Calculating error estimate."); 
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double> *>(&coarse_rho,&coarse_rho_v_x,&coarse_rho_v_y,&coarse_rho_e), Hermes::vector<Solution<double> *>(&new_rho, &new_rho_v_x, &new_rho_v_y, &new_e)) * 100;

      // Report results.
		info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", ndof, ref_ndof , err_est_rel_total);

      // If err_est too large, adapt the meshes.
      if ((err_est_rel_total < ERR_STOP)||(as>ADAPSTEP_MAX))
        done = true;
      else 
      {
    info("Adapting the coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector,&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
 
			 // Visualize the solution.
			if(HERMES_VISUALIZATION){
				s1.show(&new_rho);
  		}
			if(HERMES_VISUALIZATION_ORDER)
			{
				mview_r.show(&space_rho);
				mview_vx.show(&space_rho_v_x);
				mview_vy.show(&space_rho_v_y);
				mview_e.show(&space_e);
			}	
  		
			if(done==true)
			{
				prev_rho.copy(&new_rho); 
				prev_rho_v_x.copy(&new_rho_v_x); 
				prev_rho_v_y.copy(&new_rho_v_y); 
				prev_e.copy(&new_e); 
				

				// Output solution in VTK format.
				if((VTK_VISUALIZATION)&&(ts  % VTK_FREQ == 0)) 
				{
					sprintf(filename, "pressure-%i.vtk", ts);
					lin_p.save_solution_vtk(&pressure, filename, "Pressure", true);  
					sprintf(filename, "density-%i.vtk", ts);
					lin_rho.save_solution_vtk(&prev_rho, filename, "density", true); 
					sprintf(filename, "mesh-%i.vtk", ts);
					ord_p.save_mesh_vtk((*ref_spaces)[0], filename);   
				}			

			}	





			  // Clean up.
		delete ref_spaces;
		delete new_rho.get_mesh();
		delete new_rho_v_x.get_mesh();
		delete new_rho_v_y.get_mesh();
		delete new_e.get_mesh();
		delete lowOrd;
		delete diffusion;
		low_matrix->free();
		matrix_L->free();
		delete lumped_matrix;
		delete vec_rhs;
		delete [] P_plus;
		delete [] P_minus;
		delete [] Q_plus;
		delete [] Q_minus;
		delete [] R_plus;
		delete [] R_minus;
		delete [] coeff_vec_2;
		delete [] coeff_vec; 
		delete dp_mass;
		delete dp_boundary;
		delete dp_K;
				delete diffusion_low;
		delete dp_boundary_low;
		delete dp_K_low;


		as++;

    }
    while (done == false);

			 // Visualize the solution.
			if(HERMES_VISUALIZATION){
				sprintf(title, "pressure: ts=%i",ts);
				pressure_view.set_title(title);
				s1.show(&prev_rho);
				s2.show(&vel_x);
				s3.show(&vel_y);
				pressure_view.show(&pressure);
  		}






		
		
		
	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;




}
while (current_time < T_FINAL);

     if(VTK_VISUALIZATION)
     {
			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho_end.vtk", "density", true);
			ord_p.save_mesh_vtk(&space_rho, "mesh_end.vtk");
		}


		//Cleanup
			delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete low_matrix;
			delete lowmat_rhs;


		

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

