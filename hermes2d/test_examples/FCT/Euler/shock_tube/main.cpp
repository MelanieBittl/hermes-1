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



const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 1e-3;
const double T_FINAL = 0.231;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;  

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 



//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"


     

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);


SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));	


	int dof_rho = space_rho->get_num_dofs();
	int dof_v_x = space_rho_v_x->get_num_dofs();
	int dof_v_y = space_rho_v_y->get_num_dofs();
	int dof_e = space_e->get_num_dofs();

    Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  printf("ndof: %d \n", ndof);



  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));	

    MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_e(new Solution<double>);

  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> boundary_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	

    MeshFunctionSharedPtr<double> low_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double>	low_e(new Solution<double>);





 //--------- Visualization of pressure & velocity
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));
      s1.set_min_max_range(0, 1.);
      s2.set_min_max_range(0., 1.);
      s3.set_min_max_range(0., 1.);
			pressure_view.set_min_max_range(0.,1.);

//------------


  EulerEquationsWeakForm_K  wf_K_init(KAPPA, time_step, init_rho, init_rho_v_x, init_rho_v_y, init_e);
  EulerBoundary wf_boundary_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);

  EulerEquationsWeakForm_Mass wf_mass(time_step);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, time_step, low_rho, low_rho_v_x, low_rho_v_y, low_e);
  EulerBoundary wf_boundary(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerBoundary wf_boundary_low(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, low_rho, low_rho_v_x, low_rho_v_y, low_e);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_boundary_init(&wf_boundary_init,spaces);
  DiscreteProblem<double> dp_K_init(&wf_K_init, spaces);
    
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_boundary(&wf_boundary, spaces);
  DiscreteProblem<double> dp_boundary_low(&wf_boundary_low, spaces);

  DiscreteProblem<double> dp_K(&wf_K, spaces);
  DiscreteProblem<double> dp_K_low(&wf_K_low, spaces);

  // Set up the solver, matrix, and rhs according to the solver selection.
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   
	UMFPackMatrix<double> * matrix_L_low = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS_low = new UMFPackMatrix<double> ; 


			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ndof);
				double* u_L = NULL; 
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	

    OGProjection<double> ogProjection;
		Lumped_Projection lumpedProjection;
		
    dp_mass.assemble(mass_matrix);

	//----------------------MassLumping M_L--------------------------------------------------------------------
		UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of the initial condition
	lumpedProjection.project_lumped(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec, matrix_solver);

   ogProjection.project_global(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec_2,  Hermes::vector<NormType> (HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM));

		lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

				Solution<double>::vector_to_solutions(coeff_vec, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));	



			 // Visualize the solution.
			//s1.show(prev_rho);
//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];


//Timestep loop
do
{	 

  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time);
Space<double>::assign_dofs(spaces);
 	  
 	  if(ts!=1)
 	  {
	 	  //Hermes::Mixins::Loggable::Static::info("assemble dS");
			dp_boundary.assemble(matrix_dS);
	 		//Hermes::Mixins::Loggable::Static::info("assemble K");
		  dp_K.assemble(lowmat_rhs);
		}else{
			//Hermes::Mixins::Loggable::Static::info("assemble dS");
			dp_boundary_init.assemble(matrix_dS);
	 		//Hermes::Mixins::Loggable::Static::info("assemble K");
		  dp_K_init.assemble(lowmat_rhs);		
		}


					//------------------------artificial DIFFUSION D---------------------------------------		
					  //	Hermes::Mixins::Loggable::Static::info("artificial Diffusion");
			UMFPackMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);

//	Hermes::Mixins::Loggable::Static::info("K+D");
			lowmat_rhs->add_matrix(diffusion); //L(U)=K+D
	//Hermes::Mixins::Loggable::Static::info("L+S");
			lowmat_rhs->add_matrix(matrix_dS); //L(U)+dS(U) 

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)

			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)



	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 


	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
				//		  	Hermes::Mixins::Loggable::Static::info("Solution low order ");
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
  try
  {
   lowOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(low_rho,low_rho_v_x,low_rho_v_y,low_e));

	
			dp_boundary_low.assemble(matrix_dS_low);	
    	dp_K_low.assemble(matrix_L_low);
			UMFPackMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_matrix(diffusion_low); //L(U)
			matrix_L_low->add_matrix(matrix_dS_low); //L(U)+dS(U) 

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		//Hermes::Mixins::Loggable::Static::info("antidiffusive fluxes ");
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

		for(int i=0; i<ndof;i++)
						 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);					

				Solution<double>::vector_to_solutions(coeff_vec, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));	



// PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
//VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
 // VelocityFilter vel_y(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 2);

			 // Visualize the solution.
			 sprintf(title, " ts=%i",ts);
			// pressure_view.set_title(title);
			 s1.set_title(title);
			s1.show(prev_rho);
			//s2.show(&vel_x);
			//s3.show(&vel_y);
  		//pressure_view.show(&pressure);
  		

//	View::wait(HERMES_WAIT_KEYPRESS);



 /* // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) 
    {
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) 
      {
        pressure.reinit();
        Mach_number.reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin_pressure.save_solution_vtk(&pressure, filename, "Pressure", true);

      }
    }*/




		delete lowOrd;
		delete diffusion;
		delete diffusion_low;
		low_matrix->free();


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);

/*  PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 2);
        Linearizer lin_p;
			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
        Linearizer lin_v_x;
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
        Linearizer lin_v_y;
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_rho, "rho_end.vtk", "density", true);
*/

		//Cleanup
		delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete lumped_matrix;
			delete low_matrix;
			delete lowmat_rhs;
			delete vec_rhs;
			  // Clean up.
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete[] coeff_vec_2;
			delete [] coeff_vec; 



  // Wait for the view to be closed.
  View::wait();
  return 0;
}

