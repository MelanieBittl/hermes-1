#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "interface.h"
#include "lumped_projection.h"
//#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;



const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
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




     

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain2.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));	


SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));*/


SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(mesh, P_INIT));


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
  EulerBoundary wf_boundary(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);


EulerInterface wf_DG_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e);
EulerInterface wf_DG(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_boundary_init(&wf_boundary_init,spaces);
  DiscreteProblem<double> dp_K_init(&wf_K_init, spaces);
    
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_boundary(&wf_boundary, spaces);
  DiscreteProblem<double> dp_K(&wf_K, spaces);

  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);


  // Set up the solver, matrix, and rhs according to the solver selection.
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;
	CSCMatrix<double> * matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * mat_rhs = new CSCMatrix<double> ; 
	CSCMatrix<double> * matrix_dS = new CSCMatrix<double> ; 
	CSCMatrix<double> * matrix_DG = new CSCMatrix<double> ; 



			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
				double* u_L = NULL; 
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	

    OGProjection<double> ogProjection;
	Lumped_Projection lumpedProjection;



//Projection of the initial condition
 ogProjection.project_global(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec, Hermes::vector<NormType> (HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM));
//	lumpedProjection.project_lumped(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec, matrix_solver);

/*
SpaceSharedPtr<double> space_rho_proj(new H1Space<double>(mesh, 1));	
SpaceSharedPtr<double> space_rho_v_x_proj(new H1Space<double>(mesh, 1));	
SpaceSharedPtr<double> space_rho_v_y_proj(new H1Space<double>(mesh, 1));	
SpaceSharedPtr<double> space_e_proj(new H1Space<double>(mesh, 1));
    Hermes::vector<SpaceSharedPtr<double> > spaces_proj(space_rho_proj, space_rho_v_x_proj, space_rho_v_y_proj, space_e_proj);
Space<double>::assign_dofs(spaces_proj);
   int ndof_proj = Space<double>::get_num_dofs(spaces_proj);
			double* coeff_vec_proj = new double[ndof];
	lumpedProjection.project_lumped(spaces_proj,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec_proj, matrix_solver);

				Solution<double>::vector_to_solutions(coeff_vec_proj, spaces_proj, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));	

ogProjection.project_global(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), coeff_vec, Hermes::vector<NormType> (HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM));*/

				Solution<double>::vector_to_solutions(coeff_vec, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));	

			 // Visualize the solution.
			s1.show(prev_rho);
View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
 PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 2);


Space<double>::assign_dofs(spaces);
    dp_mass.assemble(mass_matrix);
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
		  dp_K.assemble(mat_rhs);
		  dp_DG.assemble(matrix_DG);	
		}else{
			//Hermes::Mixins::Loggable::Static::info("assemble dS");
			dp_boundary_init.assemble(matrix_dS);
	 		//Hermes::Mixins::Loggable::Static::info("assemble K");
		  dp_K_init.assemble(mat_rhs);	
		  dp_DG_init.assemble(matrix_DG);	
		}


					//------------------------artificial DIFFUSION D---------------------------------------		

			mat_rhs->add_sparse_matrix(matrix_dS); //L(U)+dS(U) 
			mat_rhs->add_sparse_matrix(matrix_DG);

			matrix->create(mat_rhs->get_size(),mat_rhs->get_nnz(), mat_rhs->get_Ap(), mat_rhs->get_Ai(),mat_rhs->get_Ax());
			matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			matrix->add_sparse_matrix(mass_matrix); 				//M/t - theta L(U)

			mat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			mat_rhs->add_sparse_matrix(mass_matrix);  //M_L/t+(1-theta)L(U)



	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			mat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 


	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
				//		  	Hermes::Mixins::Loggable::Static::info("Solution low order ");
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
  try
  {
   solver->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_L = solver->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));

	
for(int i=0; i<ndof; i++) coeff_vec[i] = u_L[i];





			 // Visualize the solution.
			 sprintf(title, " ts=%i",ts);
			// pressure_view.set_title(title);
			 s1.set_title(title);
			s1.show(prev_rho);
		//	s2.show(&vel_x);
		//	s3.show(&vel_y);
  	//	pressure_view.show(&pressure);
  		

	//View::wait(HERMES_WAIT_KEYPRESS);



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




		delete solver;
		matrix->free();


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
			delete matrix_dS;	
			delete matrix;
			delete mat_rhs;
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

