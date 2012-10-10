#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "lumped_projection.h"
#include <list>

#define _USE_MATH_DEFINES
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;



const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 0.2;                       // Time interval length. 

const double theta = 0.5;

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
const unsigned int EVERY_NTH_STEP = 100;

int main(int argc, char* argv[])
{

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();


   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

	 double diag;			
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h = diag;
	info("h_min=%f", h);




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
  CustomInitialCondition_v_x prev_rho_v_x(&mesh);
  CustomInitialCondition_v_y prev_rho_v_y(&mesh);
CustomInitialCondition_e prev_e(&mesh, KAPPA);

 Solution<double> high_rho,high_rho_v_x,high_rho_v_y,high_rho_e;

//InitialCondition auch fuer low setzen wegen Filtern!
	CustomInitialCondition_rho low_rho(&mesh);
  CustomInitialCondition_v_x low_rho_v_x(&mesh);
  CustomInitialCondition_v_y low_rho_v_y(&mesh);
CustomInitialCondition_e low_rho_e(&mesh, KAPPA);

	double current_time = 0.0; 
	CustomBoundaryCondition_rho boundary_rho(&mesh, current_time);
  CustomBoundaryCondition_v_x boundary_v_x(&mesh,  current_time);
  CustomBoundaryCondition_v_y boundary_v_y(&mesh, current_time);
		CustomBoundaryCondition_e boundary_e(&mesh, current_time, KAPPA);



 //--------- Filters for visualization of pressure & velocity
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 2);
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));
      //  s1.set_min_max_range(1.4, 8.);
	/*		s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);
/*	PressureFilter pressure_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), KAPPA);
	VelocityFilter vel_x_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 1);
	VelocityFilter vel_y_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 2);
		ScalarView s1_n("low_rho", new WinGeom(0, 0, 600, 300));
		ScalarView s2_n("low_rho_v_x", new WinGeom(700, 0, 600, 300));
		ScalarView s3_n("low_rho_v_y", new WinGeom(0, 400, 600, 300));
		ScalarView s4_n("low_pressure", new WinGeom(700, 400, 600, 300));
     s1_n.set_min_max_range(0., 1.);
      s2_n.set_min_max_range(0., 1.);
      s3_n.set_min_max_range(0., 1.);
      s4_n.set_min_max_range(0., 1.);*/
      
      
			Linearizer lin_p;
			Linearizer lin_v_x;
			Linearizer lin_v_y;
			Linearizer lin_rho;

		/*	lin_p.save_solution_vtk(&pressure, "p_init.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "vx_init.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_init.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho_init.vtk", "density", true);*/
      
      
//------------

  EulerEquationsWeakForm_Mass wf_mass(time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, time_step, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);


  EulerBoundary wf_boundary(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerBoundary wf_boundary_low(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  DiscreteProblem<double> dp_boundary(&wf_boundary, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_boundary_low(&wf_boundary_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  DiscreteProblem<double> dp_K(&wf_K, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_K_low(&wf_K_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));


  // Set up the solver, matrix, and rhs according to the solver selection.
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   
	UMFPackMatrix<double> * matrix_L_low = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 

	UMFPackMatrix<double> * matrix_dS = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS_low = new UMFPackMatrix<double> ; 

	UMFPackVector<double> * surf_rhs= new UMFPackVector<double> (ndof);;

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ndof);
				double* u_L = NULL; 
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	




    dp_mass.assemble(mass_matrix);
	//----------------------MassLumping M_L--------------------------------------------------------------------
		UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of the initial condition
			Lumped_Projection::project_lumped(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec, matrix_solver);

			OGProjection<double>::project_global(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec_2, matrix_solver, HERMES_L2_NORM);

			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

/*			Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));	

			s1_n.show(&low_rho);
			s2_n.show(&vel_x_low);
			s3_n.show(&vel_y_low);
			s4_n.show(&pressure_low);

		/*Solution<double>::vector_to_solutions(coeff_vec_2, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&high_rho,&high_rho_v_x,&high_rho_v_y,&high_rho_e));	
			s1_n.show(&high_rho);
			s2_n.show(&high_rho_v_x);
			s3_n.show(&high_rho_v_y);
			s4_n.show(&high_rho_e);*/

	//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:

	int ts = 1;
	char title[100];




//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

// Time measurement.
  cpu_time.tick(HERMES_SKIP);

		dp_boundary.assemble(matrix_dS,surf_rhs);
    dp_K.assemble(lowmat_rhs);

					//------------------------artificial DIFFUSION D---------------------------------------		
			UMFPackMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);
			
			lowmat_rhs->add_matrix(diffusion); //L(U)
			lowmat_rhs->add_matrix(matrix_dS); //L(U)+dS(U) 

//matrix_L = K+D+dS
			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)
			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)


	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			
			for(int i =0;i<ndof;i++) coeff_vec_2[i] += surf_rhs->get(i)*time_step;
						vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 


	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));
			  }else error("Matrix solver failed.\n");



		dp_boundary_low.assemble(matrix_dS_low);	
    dp_K_low.assemble(matrix_L_low);
			UMFPackMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_matrix(diffusion_low); //L(U)	
			matrix_L_low->add_matrix(matrix_dS_low); //L(U)+dS(U) 


		//---------------------------------------antidiffusive fluxes-----------------------------------	
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low,surf_rhs,  u_L,coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);
				for(int i=0; i<ndof;i++){
								 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		
				}

				Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));	



			 // Visualize the solution.
	 sprintf(title, "pressure: ts=%i",ts);
			 pressure_view.set_title(title);
			s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);


	//View::wait(HERMES_WAIT_KEYPRESS);



  // Visualization.
   if((ts - 1) % EVERY_NTH_STEP == 0) 
    {
      // Output solution in VTK format.

        char filename[40];
       /* sprintf(filename, "pressure-%i.vtk", ts - 1);
        lin_p.save_solution_vtk(&pressure, filename, "Pressure", true);
        			        sprintf(filename, "v_x-%i.vtk", ts - 1);
			lin_v_x.save_solution_vtk(&vel_x, filename,  "velocity", true); */
        sprintf(filename, "rho-%i.vtk", ts - 1);
			lin_rho.save_solution_vtk(&prev_rho, filename,  "density", true);

      
    }




		delete lowOrd;
		delete diffusion;
		delete diffusion_low;
		low_matrix->free();


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

boundary_rho.set_time(current_time);
boundary_v_x.set_time(current_time);
boundary_v_y.set_time(current_time);
boundary_e.set_time(current_time);

}
while (current_time < T_FINAL);

			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho_end.vtk", "density", true);


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

