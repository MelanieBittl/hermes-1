#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"

#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 


const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
                    
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length. 
const double EPS = 1e-6;

const int NDOF_STOP = 20000;   

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

//FCT & Adaptivity....
#include "fct.cpp"

#include "h_adapt.cpp"
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "p1_list.cpp"
#include "reg_estimator.cpp"




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


  // Initialize boundary conditions.
  DefaultEssentialBCConst<double>  bc_essential(BDY_IN, 0.0);
  EssentialBCs<double>  bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

 // Initialize solution of lower & higher order
  Solution<double>  low_sln, u_new;
		Solution<double> R_h_1, R_h_2;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  

	//CustomInitialCondition p1_sln(&mesh);  	
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, &u_prev_time);
	CustomWeakFormConvection  convection(&u_prev_time);


  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(&u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
	//mview.show(&space);


  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
lin.save_solution_vtk(&u_prev_time, "init_uniform.vtk", "u", mode_3D);

		  // Initialize
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 
	double* u_L = NULL; 
	double* u_new_double = NULL; 



// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;  //p-adapt-schritte

	//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
	int ref_ndof;


	H1Space<double>* ref_space = new H1Space<double>(&mesh, &bcs, P_INIT);	

	DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, ref_space);
	DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, ref_space);
	AsmList<double>*  al = new AsmList<double>;	

	AsmList<double>*  dof_list= new AsmList<double>;
			ref_ndof = ref_space->get_num_dofs();


			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];
			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];	
			double* Q_plus_old = new double[ref_ndof]; double* Q_minus_old = new double[ref_ndof];	
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];	
				for(int i=0; i<ref_ndof;i++){ coeff_vec[i]=0.0;	coeff_vec_2[i]=0.0;	Q_plus_old[i]=0.;	Q_minus_old[i]=0.;}

			p1_list_fast(ref_space, dof_list, al,P_plus,P_minus);


				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			  // Set up the matrix, and rhs according to the solver selection.=>For Masslumping
		//	info("Masslumping");

			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			  // Set up the solver, matrix, and rhs according to the solver selection.=>artificial Diffusion
			//info("artificialDiffusion");

			dp_convection->assemble(conv_matrix, NULL,true);
			UMFPackMatrix<double> * diffusion = artificialDiffusion(conv_matrix);

			//--------------------------------------------------------------------------------------------

			lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			lowmat_rhs->add_matrix(diffusion); 
			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			//(-theta)(K+D)
			if(theta==0) low_matrix->zero();
			else	low_matrix->multiply_with_Scalar(-theta);
			//(1-theta)(K+D)
			if(theta ==1) lowmat_rhs->zero();
			else lowmat_rhs->multiply_with_Scalar((1.0-theta));

			//M_L/tau - theta(D+K)
			low_matrix->add_matrix(lumped_matrix);  
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs->add_matrix(lumped_matrix);	

			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C

		int* smooth_fct_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dx_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dy_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[ref_space->get_num_dofs()];
			double* flux_double = new double[ref_ndof]; 
			double* lumped_double = new double[ref_ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ref_ndof);


		//Initialisierung von Q_plus_old,Q_minus_old
	/*for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
			if(mass_matrix->get(i,j)!=0.0){
						f = lumped_matrix->get(i,i)*(u_prev_time.get_pt_value(P_plus[j], P_minus[j])- u_prev_time.get_pt_value(P_plus[i], P_minus[i])); 
						if(f>Q_plus_old[i]) Q_plus_old[i] = f;				
						if(f<Q_minus_old[i]) Q_minus_old[i] = f;			
						f= lumped_matrix->get(j,j)*(u_prev_time.get_pt_value(P_plus[i], P_minus[i])- u_prev_time.get_pt_value(P_plus[j], P_minus[j])); 
						if(f>Q_plus_old[j]) Q_plus_old[j] = f;	
						if(f<Q_minus_old[j]) Q_minus_old[j] = f;
					}
				}
			}	*/



			// Project the initial condition on the FE space->coeff_vec	
			//info("projection");
			Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
			Solution<double> ::vector_to_solution(coeff_vec, ref_space, &low_sln);
		smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
			OGProjection<double>::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old,  R_plus, R_minus,smooth_dof);

			Solution<double> ::vector_to_solution(coeff_vec, ref_space, &u_new);

			/*sprintf(title, "proj. Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(&u_new);
			sprintf(title, "proj. lumpedLoesung, ps=%i, ts=%i", ps,ts);
			Lowview.show(&low_sln);
			mview.show(ref_space);		
			View::wait(HERMES_WAIT_KEYPRESS);*/
	





//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 
			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);

  // Time measurement.
 // cpu_time.tick(HERMES_SKIP);

	//-------------------------solution of lower order------------	
//info("Linear Solver lowOrd");
			  // Solve the linear system and if successful, obtain the solution. M_L/tau-theta(D+K) u^n+1=  M_L/tau+ (1-theta)(K+D) u^n
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_L, ref_space, &low_sln);	
			  }else error ("Matrix solver failed.\n");

  // CPU time needed by UMFpack
  //double time1 = cpu_time.tick().last();
  // Time measurement.
  //cpu_time.tick(HERMES_SKIP);

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
			//	info("assemble fluxes");	
			 antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,vec_rhs, u_L,coeff_vec, flux_double, 
									P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus,smooth_dof);


				for(int i=0; i<ref_ndof;i++){
								 coeff_vec[i] = u_L[i]+ (flux_double[i]*time_step)/lumped_matrix->get(i,i);
									Q_plus_old[i]=0.;Q_minus_old[i]=0.;
				}
				Solution<double>::vector_to_solution(coeff_vec, ref_space, &u_new);	

  // Measure the projection time.
 // double flux_time = cpu_time.tick().last();

//info("CPU_Time:  %g -------------- %g", time1, flux_time);

		//Initialisierung von Q_plus_old,Q_minus_old
	/*for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
			if(mass_matrix->get(i,j)!=0.0){
						f = lumped_matrix->get(i,i)*(coeff_vec[j]- coeff_vec[i])/time_step; 
						if(f>Q_plus_old[i]) Q_plus_old[i] = f;				
						if(f<Q_minus_old[i]) Q_minus_old[i] = f;			
						f= lumped_matrix->get(j,j)*(coeff_vec[i]- coeff_vec[j])/time_step; 
						if(f>Q_plus_old[j]) Q_plus_old[j] = f;	
						if(f<Q_minus_old[j]) Q_minus_old[j] = f;
					}
				}
			}	*/

			 // Visualize the solution.
			  sprintf(title, "korrigierte Loesung: ps=%i, ts=%i", ps,ts);
			  sview.set_title(title);
		  sview.show(&u_new);
				//mview.show(ref_space);
	//View::wait(HERMES_WAIT_KEYPRESS);

			

			delete lowOrd;
			//delete newSol; 	
		

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);

lin.save_solution_vtk(&u_new, "end_uniform_smooth.vtk", "solution", mode_3D);

		delete dp_convection;
		delete dp_mass; 
		delete ref_space; 

	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;

delete al;


	  // Clean up.
			delete[] flux_double; 
			delete vec_rhs;
			delete[] lumped_double;


			delete dof_list;
		

		delete [] smooth_fct_in_elem;
		delete [] smooth_dx_in_elem;
		delete [] smooth_dy_in_elem;
		delete [] smooth_elem;
		delete [] smooth_dof;


			
			  // Clean up.
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] Q_plus_old;
			delete [] Q_minus_old;
			delete [] R_plus;
			delete [] R_minus;
			delete lumped_matrix; 
			delete diffusion;
			delete[] coeff_vec_2;
			if(coeff_vec!=NULL){ delete [] coeff_vec; 	coeff_vec = NULL;}
		  low_matrix->free();
	 		lowmat_rhs->free();


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

