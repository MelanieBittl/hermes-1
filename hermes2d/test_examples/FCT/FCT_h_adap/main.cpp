#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
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
const int P_MAX = 2; 
const double h_max = 0.1;                       
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.

const double P_ADAP_TOL_EX = 0.2;   
const int 	 P_ADAP_MAX_ITER = 2;
 
const double EPS = 1e-6;

const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 20;                  // Maximum allowed number of Newton iterations.

const int NDOF_STOP = 20000;   

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

//FCT & p-Adaptivity
//#include "fct.cpp"
#include "fct_simple.cpp"
#include "h_adapt.cpp"
//#include "mass_lumping.cpp"
//#include "artificial_diffusion.cpp"
//#include "p1_list.cpp"
//#include "reg_estimator.cpp"
#include "z_z.cpp"


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
	H1Space<double>* ref_space = new H1Space<double>(&mesh, &bcs, P_INIT);	

 // Initialize solution of lower & higher order
  Solution<double>  low_sln, u_new,high_sln, u_prev;
		Solution<double> R_h_1, R_h_2;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);   	
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, &u_prev_time);
	CustomWeakFormConvection  convection(&u_prev_time);

	ConvectionForm  wf(time_step, &u_prev_time);


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
lin.save_solution_vtk(&u_prev_time, "init_h_adap_neu.vtk", "u", mode_3D);

		  // Initialize
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 

		UMFPackMatrix<double> * high_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * high_rhs = new UMFPackMatrix<double> ; 

	double* u_L = NULL; 
	double* u_new_double = NULL; 
	double* u_H =NULL;



// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;  //p-adapt-schritte



	 double diag;
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8; 
	info("h_min=%f", h_min);


		int ref_ndof;

//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 
	
	mesh.copy(&basemesh);
	ref_space->set_mesh(&mesh);	ref_space->set_uniform_order(1); ref_space->assign_dofs(); 
	ps=1; 

//Adaptivity loop
	do
	{		
			DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, ref_space);
			DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, ref_space);
			HPAdapt * adapting = new HPAdapt(ref_space, HERMES_L2_NORM);

			ref_ndof = ref_space->get_num_dofs();

			info(" adap- step %d, timestep %d,ndof = %d ", ps, ts, ref_ndof); 

			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];
			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];	
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];	
	
			double* coeff_vec_3= new double[ref_ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ref_ndof);
			UMFPackVector<double> * vec_high_rhs = new UMFPackVector<double> (ref_ndof);
				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<double> * lumped_matrix = massLumping_simple(mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection->assemble(conv_matrix, NULL,true);
			UMFPackMatrix<double> * diffusion = artificialDiffusion_simple(conv_matrix);
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

			
		high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_matrix->multiply_with_Scalar(-theta);
			high_matrix->add_matrix(mass_matrix);  
			high_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_matrix(mass_matrix); 

			lumped_matrix->multiply_with_Scalar(time_step);  // M_L

	//--------Project the initial condition on the FE space->coeff_vec	---------------

			if(ts==1) {
				Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec,  matrix_solver, lumped_matrix);
							OGProjection<double>::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			}else{
				 Lumped_Projection::project_lumped(ref_space, &u_prev, coeff_vec,  matrix_solver, lumped_matrix);
			OGProjection<double>::project_global(ref_space,&u_prev, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			}
			lumped_flux_limiter_simple(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);

		/*	Solution<double>::vector_to_solution(coeff_vec, ref_space, &u_new);

			sprintf(title, "proj. Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(&u_new);*/

		//	sprintf(title, "proj. lumpedLoesung, ps=%i, ts=%i", ps,ts);

			//mview.show(ref_space);		
			//View::wait(HERMES_WAIT_KEYPRESS);


	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_3); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_3);

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L (u^L/tau)=  M_L/tau+ (1-theta)(K+D) u^n

			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (lumped_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
			for(int i=0;i<ref_ndof;i++) u_L[i]=u_L[i]*time_step;
				Solution<double> ::vector_to_solution(u_L, ref_space, &low_sln);	
			  }else error ("Matrix solver failed.\n");


			/*	sprintf(title, "low sln Time %3.2f,timestep %i,ps=%i,", current_time,ts,ps);
				 sview.set_title(title);
					sview.show(&low_sln);
*/

		//-------------solution of higher order

			high_rhs->multiply_with_vector(coeff_vec, coeff_vec_3); 
			vec_high_rhs->zero(); vec_high_rhs->add_vector(coeff_vec_3);
			UMFPackLinearSolver<double> * highOrd = new UMFPackLinearSolver<double> (high_matrix,vec_high_rhs);	
			if(highOrd->solve()){ 
				u_H = highOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_H, ref_space, &high_sln);	
			  }else error ("Matrix solver failed.\n");

   	/*sprintf(title, "high-sln adap-step %i", ps);
			 sview.set_title(title);
			 sview.show(&high_sln);*/

		//---------------------------------------antidiffusive fluxes-----------------------------------	
antidiffusiveFlux_simple(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, u_L,coeff_vec, coeff_vec_3, P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);

			vec_rhs->add_vector(coeff_vec_3);
			UMFPackLinearSolver<double> * newSol = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(newSol->solve()){ 
				u_new_double = newSol->get_sln_vector();  
				Solution<double>::vector_to_solution(u_new_double, ref_space, &u_new);	
			}else error ("Matrix solver failed.\n");	 



			 // Visualize the solution.			 
			if(ps==2){	
				sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i,ps=%i,", current_time,ts,ps);
				 sview.set_title(title);
					sview.show(&u_new);
			}


if(ps==2) 
		u_prev.copy(&u_new);

	if(ps==1){	
	changed = h_p_adap(ref_space,&u_new,&R_h_1,&R_h_2,adapting,h_min,h_max,ps);	
	
			/*sprintf(title, "nach changed Mesh, ps=%i, ts=%i", ps,ts);
			mview.set_title(title);
				mview.show(ref_space);*/

			}
			


			ps++;	
	  // Clean up.
		delete[] coeff_vec_3; 
			delete vec_rhs;

			delete vec_high_rhs;
			high_matrix->free();
			high_rhs->free();

		delete newSol;
		delete highOrd;
		delete lowOrd;


			  // Clean up.
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete lumped_matrix; 
			delete diffusion;
			delete[] coeff_vec_2;
			 delete [] coeff_vec; 
		  low_matrix->free();
	 		lowmat_rhs->free();

		delete dp_convection;
		delete dp_mass; 
		delete adapting;



	}while(ps<3);		

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
 	ts++;

}
while (current_time < T_FINAL);

lin.save_solution_vtk(&u_prev_time, "end_h_adap_neu.vtk", "solution", mode_3D);
/*sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&u_new);*/


		delete ref_space; 

	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;
	delete high_matrix;
	delete high_rhs;





  // Wait for the view to be closed.
  View::wait();
  return 0;
}

