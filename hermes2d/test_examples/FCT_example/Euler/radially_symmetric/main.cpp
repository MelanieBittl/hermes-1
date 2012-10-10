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



const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 2e-3;
const double T_FINAL = 0.13;                       // Time interval length. 

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
const unsigned int EVERY_NTH_STEP = 2;

int main(int argc, char* argv[])
{
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
  ConstantSolution<double> prev_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, V2_EXT);
CustomInitialCondition_e prev_e(&mesh, KAPPA);


 Solution<double> high_rho,high_rho_v_x,high_rho_v_y,high_rho_e;

//InitialCondition auch fuer low setzen wegen Filtern!
	CustomInitialCondition_rho low_rho(&mesh);
  ConstantSolution<double> low_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> low_rho_v_y(&mesh, V2_EXT);
CustomInitialCondition_e low_rho_e(&mesh, KAPPA);

	CustomInitialCondition_rho boundary_rho(&mesh);
  ConstantSolution<double> boundary_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> boundary_v_y(&mesh, V2_EXT);
CustomInitialCondition_e boundary_e(&mesh,KAPPA);

        Linearizer lin_p;
        Linearizer lin_v_x; 
        Linearizer lin_v_y;
        Linearizer lin_rho;

 //--------- Filters for visualization of pressure & velocity
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 2);
  RadiusVelocityFilter vel_r(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));
  ScalarView radius_view("Radius", new WinGeom(700, 700, 600, 300));
      s1.set_min_max_range(0, 2.5);
			pressure_view.set_min_max_range(0.,4.);
		/*	s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);*/
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
      
			lin_p.save_solution_vtk(&pressure, "pressure-0.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "v_x-0.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy-0.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho-0.vtk", "density", true);
      
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
	double current_time = 0.0; 
	int ts = 1;
	char title[100];




//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 


		dp_boundary.assemble(matrix_dS);
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
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);
				for(int i=0; i<ndof;i++){
								 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		
				}

				Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));	



			 // Visualize the solution.
	/*	  sprintf(title, "pressure: ts=%i",ts);
			 pressure_view.set_title(title);
			s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);
			radius_view.show(&vel_r);*/


//	View::wait(HERMES_WAIT_KEYPRESS);


// Visualization.
  if((ts - 1) % EVERY_NTH_STEP == 0) 
    {
      // Output solution in VTK format.
       // pressure.reinit();
        char filename[40];
        sprintf(filename, "pressure-%i.vtk", ts );
        lin_p.save_solution_vtk(&pressure, filename, "Pressure", true);
        sprintf(filename, "rho-%i.vtk", ts );
			lin_rho.save_solution_vtk(&prev_rho, filename,  "density", true);
			        sprintf(filename, "v_x-%i.vtk", ts );
			lin_v_x.save_solution_vtk(&vel_x, filename,  "density", true);
      
    }




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

			lin_p.save_solution_vtk(&pressure, "pressure-130.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "v_x-130.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho-130.vtk", "density", true);


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

