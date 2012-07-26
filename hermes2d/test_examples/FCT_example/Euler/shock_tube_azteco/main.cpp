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

MatrixSolverType matrix_solver = SOLVER_AZTECOO; 

//MatrixSolverType matrix_solver = SOLVER_UMFPACK; 



//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"


     

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

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
		/*	s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);*/
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

  EulerEquationsWeakForm_Mass wf_mass(time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, time_step, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);
  EulerEquationsWeakForm_Surf  wf_surf(KAPPA, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  EulerBoundary wf_boundary(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerBoundary wf_boundary_low(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  DiscreteProblem<double> dp_boundary(&wf_boundary, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_boundary_low(&wf_boundary_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  DiscreteProblem<double> dp_K(&wf_K, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
    DiscreteProblem<double> dp_K2(&wf_K, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_K_low(&wf_K_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));


  // Set up the solver, matrix, and rhs according to the solver selection.
	EpetraMatrix<double> * mass_matrix = new EpetraMatrix<double> ;   
	EpetraMatrix<double> * matrix_L_low = new EpetraMatrix<double> ; 
	EpetraMatrix<double> * low_matrix = new EpetraMatrix<double> ;  
	EpetraMatrix<double> * lowmat_rhs = new EpetraMatrix<double> ; 

	EpetraMatrix<double> * matrix_dS = new EpetraMatrix<double> ; 
	EpetraMatrix<double> * matrix_dS_low = new EpetraMatrix<double> ; 



			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			EpetraVector<double> * vec_rhs = new EpetraVector<double>;
 				vec_rhs->alloc(ndof);

				double* u_L = NULL; 
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	


    dp_mass.assemble(mass_matrix);

	//----------------------MassLumping M_L--------------------------------------------------------------------
		EpetraMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of the initial condition
//Lumped Projection noch mit UMFPack!!
			Lumped_Projection::project_lumped(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec);

			OGProjection<double>::project_global(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec_2, matrix_solver, HERMES_L2_NORM);


		lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

	/*Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));	

			s1.show(&low_rho);
			s2.show(&vel_x_low);
			s3.show(&vel_y_low);
			s4.show(&pressure_low);

	/*	Solution<double>::vector_to_solutions(coeff_vec_2, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&high_rho,&high_rho_v_x,&high_rho_v_y,&high_rho_e));	
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
// Time measurement.
  cpu_time.tick(HERMES_SKIP);
  
		dp_boundary.assemble(matrix_dS);
    dp_K.assemble(lowmat_rhs);

  // CPU time 
  //double time1 = cpu_time.tick().last();
   //Time measurement.
  //cpu_time.tick(HERMES_SKIP);

					//------------------------artificial DIFFUSION D---------------------------------------		
			EpetraMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);			

							
						lowmat_rhs->add_matrix(diffusion);
						lowmat_rhs->add_matrix(matrix_dS);
			
						  // CPU time 
				double time2 = cpu_time.tick().last();
				 //Time measurement.
				cpu_time.tick(HERMES_SKIP);
			
			    dp_K2.assemble(low_matrix);
			    low_matrix->add_matrix(diffusion);
					low_matrix->add_matrix(matrix_dS);			

//matrix_L = K+D+dS
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)
			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)

lowmat_rhs->finish();
low_matrix->finish();
	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 

  // CPU time 
  //double time3 = cpu_time.tick().last();
   //Time measurement.
 // cpu_time.tick(HERMES_SKIP);

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			AztecOOSolver<double>* lowOrd = new AztecOOSolver<double> (low_matrix,vec_rhs);	
			
			      /// @param[in] solver - name of the solver [ gmres | cg | cgs | tfqmr | bicgstab ]
			lowOrd->set_solver("bicgstab");
      lowOrd->set_tolerance(1e-8);
     // lowOrd->set_max_iters(int iters);

      /// @param[in] name - name of the preconditioner [ none | jacobi | neumann | least-squares ]
     // lowOrd->set_precond("jacobi");


			if(lowOrd->solve(coeff_vec)){ 
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));
			  }else error("Matrix solver failed.\n");



			/*s1_n.show(&low_rho);
			s2_n.show(&vel_x_low);
			s3_n.show(&vel_y_low);
			s4_n.show(&pressure_low);	
	
	View::wait(HERMES_WAIT_KEYPRESS);	*/


	dp_boundary_low.assemble(matrix_dS_low);	
    dp_K_low.assemble(matrix_L_low);
			EpetraMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_matrix(diffusion_low); //L(U)	
			matrix_L_low->add_matrix( matrix_dS_low); //L(U)+dS(U) 


		//---------------------------------------antidiffusive fluxes-----------------------------------	
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);
				for(int i=0; i<ndof;i++){
								 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		
				}

				Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));	

  // Measure the solver time.
  double time4 = cpu_time.tick().last();
//info("CPU_Time:  %g -------%g------- %g----------%g", time1,time2, time3, time4);
info("CPU_Time:  %g ", time4);

			 // Visualize the solution.
		  sprintf(title, "pressure: ts=%i",ts);
			 pressure_view.set_title(title);
			s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
  		pressure_view.show(&pressure);


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




		delete lowOrd;
		delete diffusion;
		delete diffusion_low;
	


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);
        Linearizer lin_p;
			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
        Linearizer lin_v_x;
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
        Linearizer lin_v_y;
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(&prev_rho, "rho_end.vtk", "density", true);


		//Cleanup
			delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete lumped_matrix;
			delete lowmat_rhs;
			delete vec_rhs;
			delete low_matrix;
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

