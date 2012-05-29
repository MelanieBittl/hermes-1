#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "lumped_projection.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;



const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 1e-3;
const double T_FINAL = 0.3;                       // Time interval length. 

const double theta = 1.;

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_Solid = "solid";


//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"


// Equation parameters.      
// Inlet x-velocity 
const double V1_EXT = 0.0;       
// Inlet y-velocity 
const double V2_EXT = 0.0;        
// Kappa. => for air
const double KAPPA = 1.4;      

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

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


  // Initialize boundary conditions.
 // DefaultEssentialBCConst<double>  bc_essential(BDY, 0.0);
//  EssentialBCs<double>  bcs(&bc_essential);

/*    EssentialBCs<double> bcs_rho;
  bcs_rho.add_boundary_condition(new CustomBC_rho(BDY_Solid,KAPPA));
    EssentialBCs<double> bcs_e;
  bcs_e.add_boundary_condition(new CustomBC_e(BDY_Solid,KAPPA));

 DefaultEssentialBCConst<double>  bc_essential(BDY_Solid, 0.0);
EssentialBCs<double>  bcs_v_y(&bc_essential);
EssentialBCs<double>  bcs_v_x(&bc_essential);

  /*  EssentialBCs<double> bcs_v_x;
  bcs_v_x.add_boundary_condition(new CustomBC_v_x(BDY_Solid,KAPPA));
    EssentialBCs<double> bcs_v_y;
  bcs_v_y.add_boundary_condition(new CustomBC_v_y(BDY_Solid,KAPPA));*/


  // Initialize boundary condition types and spaces with default shapesets.
 /*H1Space<double> space_rho(&mesh,&bcs_rho,  P_INIT);
  H1Space<double> space_rho_v_x(&mesh,&bcs_v_x,  P_INIT);
  H1Space<double> space_rho_v_y(&mesh,&bcs_v_y,  P_INIT);
  H1Space<double> space_e(&mesh,&bcs_e,  P_INIT);*/



 H1Space<double> space_rho(&mesh,P_INIT);
  H1Space<double> space_rho_v_x(&mesh,P_INIT);
  H1Space<double> space_rho_v_y(&mesh,P_INIT);
  H1Space<double> space_e(&mesh, P_INIT);


  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
	CustomInitialCondition_rho prev_rho(&mesh);
  ConstantSolution<double> prev_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, V2_EXT);
CustomInitialCondition_e prev_e(&mesh, KAPPA);

  Solution<double>  low_rho,low_rho_v_x,low_rho_v_y,low_rho_e;
 Solution<double> high_rho,high_rho_v_x,high_rho_v_y,high_rho_e  ;
 Solution<double> new_rho,new_rho_v_x,new_rho_v_y,new_rho_e  ;

  // Filters for visualization of pressure
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);


  ScalarView pressure_view("Pressure", new WinGeom(700, 700, 600, 300));
  ScalarView s1("prev_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("prev_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("prev_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 400, 600, 300));

  ScalarView s1_n("low_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2_n("low_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3_n("low_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4_n("low_e", new WinGeom(700, 400, 600, 300));

		/*	s1.show(&prev_rho);
			s2.show(&prev_rho_v_x);
			s3.show(&prev_rho_v_y);
			s4.show(&prev_e);
  		pressure_view.show(&pressure);*/


  EulerEquationsWeakForm_Mass wf_mass(time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_Surf  wf_surf(KAPPA, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_surf(&wf_surf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_K(&wf_K, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));


  // Set up the solver, matrix, and rhs according to the solver selection.
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   
	UMFPackMatrix<double> * matrix_K = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_L = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_S = new UMFPackMatrix<double> ; 
//	UMFPackVector<double> * rhs_surf = new UMFPackVector<double> (ndof);  

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ndof);
				double* u_L = NULL; 
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	




    dp_mass.assemble(mass_matrix);
   dp_surf.assemble(matrix_S);
    dp_K.assemble(matrix_K);

	//----------------------MassLumping M_L--------------------------------------------------------------------
		UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);
		int size = matrix_K->get_size();
					//------------------------artificial DIFFUSION D---------------------------------------		
			UMFPackMatrix<double> * diffusion = artificialDiffusion(KAPPA,&prev_rho, &prev_rho_v_x, 
    &prev_rho_v_y, &prev_e,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e);

		lumped_matrix->multiply_with_Scalar(1./time_step); //M_L/tau

			lowmat_rhs->create(matrix_K->get_size(),matrix_K->get_nnz(), matrix_K->get_Ap(), matrix_K->get_Ai(),matrix_K->get_Ax());
			lowmat_rhs->add_matrix(diffusion); //L(U)
			lowmat_rhs->add_matrix(matrix_S); //L(U)+dS(U) 
			matrix_L->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)

			lowmat_rhs->multiply_with_Scalar((1.0-theta));  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)




//Projection of the initial condition
			Lumped_Projection::project_lumped(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec, matrix_solver);
			OGProjection<double>::project_global(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);

				Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&new_rho,&new_rho_v_x,&new_rho_v_y,&new_rho_e));	

			//	Solution<double>::vector_to_solutions(coeff_vec_2, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&high_rho,&high_rho_v_x,&high_rho_v_y,&high_rho_e));	

		s1_n.show(&new_rho);
			s2_n.show(&new_rho_v_x);
			s3_n.show(&new_rho_v_y);
			s4_n.show(&new_rho_e);

			/*s1_n.show(&high_rho);
			s2_n.show(&high_rho_v_x);
			s3_n.show(&high_rho_v_y);
			s4_n.show(&high_rho_e);*/

	View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];


		lumped_matrix->multiply_with_Scalar(time_step); 

//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2);

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));
			  }else error ("Matrix solver failed.\n");


			/*s1_n.show(&low_rho);
			s2_n.show(&low_rho_v_x);
			s3_n.show(&low_rho_v_y);
			s4_n.show(&low_rho_e);*/

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		antidiffusiveFlux(mass_matrix,lumped_matrix,matrix_K,diffusion, matrix_L, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus);
				for(int i=0; i<ndof;i++){
								 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		//time_step?
				}

				Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));	



			 // Visualize the solution.
			  sprintf(title, "pressure: ts=%i",ts);
			 pressure_view.set_title(title);
			s1.show(&prev_rho);
			s2.show(&prev_rho_v_x);
			s3.show(&prev_rho_v_y);
			s4.show(&prev_e);
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

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);

//Cleanup
delete mass_matrix;
delete matrix_K;
delete matrix_L;
delete lumped_matrix;
delete diffusion;
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

