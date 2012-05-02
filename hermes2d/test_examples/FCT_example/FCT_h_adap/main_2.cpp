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
#include "fct.cpp"
#include "h_adapt.cpp"
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "p1_list.cpp"
#include "reg_estimator.cpp"
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
  H1Space<double> space(&mesh, &bcs, P_INIT);

  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

 // Initialize solution of lower & higher order
  Solution<double>  low_sln, u_new, e_h,e_h_p, p1_sln, high_sln;
		Solution<double> R_h_1, R_h_2;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  

	//CustomInitialCondition p1_sln(&mesh);  	
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, &u_prev_time);
	CustomWeakFormConvection  convection(&u_prev_time);

	 ResidualForm  residual(time_step, &u_prev_time, &p1_sln);




  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(&u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
	//mview.show(&space);


  // Output solution in VTK format.
//Linearizer lin;
//bool mode_3D = true;
//lin.save_solution_vtk(&u_prev_time, "init_hpadap_neu.vtk", "u", mode_3D);

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

	//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 

	
	std::list<int> p1_elements;
	std::list<int> neighbor;



	 double diag;
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8; 
	info("h_min=%f", h_min);



	AsmList<double>*  al = new AsmList<double>;	

	H1Space<double>* ref_space = new H1Space<double>(&mesh, &bcs, P_INIT);	
	DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, ref_space);
	DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, ref_space);
	HPAdapt * adapting = new HPAdapt(ref_space, HERMES_L2_NORM);

/*	AsmList<double>*  dof_list= new AsmList<double>;
	AsmList<double>*  dof_list_2= new AsmList<double>;
	p1_list_fast(ref_space, dof_list,dof_list_2, al);

	changed = init_adapt(ref_space, &u_prev_time,&e_h,&R_h_1,&R_h_2,&massmatrix, adapting,dof_list,dof_list_2, al, &p1_elements,&neighbor,h_min,h_max);
	delete dof_list; delete dof_list_2;*/
				mview.show(ref_space);

//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

	ps=1; 

//Adaptivity loop
	do
	{			
				AsmList<double>*  dof_list= new AsmList<double>;
				AsmList<double>*  dof_list_2= new AsmList<double>;

			//if((ps==1))
					p1_list_fast(ref_space, dof_list,dof_list_2, al);
		//	else
					 //p1_list(ref_space, dof_list,dof_list_2, al,&p1_elements,&neighbor);

			int ref_ndof = ref_space->get_num_dofs();
			dp_mass->set_spaces(ref_space);
			dp_convection->set_spaces(ref_space);

			info(" adap- step %d, timestep %d,ndof = %d ", ps, ts, ref_ndof); 

			double* coeff_vec = new double[ref_ndof];
				for(int i=0; i<ref_ndof;i++) coeff_vec[i]=0.0;	
			double* coeff_vec_2 = new double[ref_ndof];
				for(int i=0; i<ref_ndof;i++) coeff_vec_2[i]=0.0;
			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];	
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];	
				int* smooth_fct_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dx_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dy_in_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[ref_ndof];


				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			  // Set up the matrix, and rhs according to the solver selection.=>For Masslumping
		//	info("Masslumping");

			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<double> * lumped_matrix = massLumping(dof_list,dof_list_2,mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			  // Set up the solver, matrix, and rhs according to the solver selection.=>artificial Diffusion
			//info("artificialDiffusion");

			dp_convection->assemble(conv_matrix, NULL,true);
			UMFPackMatrix<double> * diffusion = artificialDiffusion(dof_list,dof_list_2,conv_matrix);

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
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
		
			// Project the initial condition on the FE space->coeff_vec	
			//info("projection");
			Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
					//	Solution<double> ::vector_to_solution(coeff_vec, ref_space, &low_sln);
		//smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
			OGProjection<double>::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(dof_list,dof_list_2,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);
			Solution<double>::vector_to_solution(coeff_vec, ref_space, &u_new);

			sprintf(title, "proj. Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(&u_new);
		//	sprintf(title, "proj. lumpedLoesung, ps=%i, ts=%i", ps,ts);
		//	Lowview.show(&low_sln);
			//mview.show(ref_space);		
			//View::wait(HERMES_WAIT_KEYPRESS);
	
			double* flux_double = new double[ref_ndof]; 
			double* lumped_double = new double[ref_ndof];
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ref_ndof);

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 
			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
			lumped_matrix->multiply_with_Scalar(1./time_step); //M_L/tau
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (lumped_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_L, ref_space, &low_sln);	
			  }else error ("Matrix solver failed.\n");
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L


		//-------------solution of higher order

	double* high_double = new double[ref_ndof];  memset(high_double, 0, ref_ndof*sizeof(double));
		UMFPackVector<double> * vec_high_rhs = new UMFPackVector<double> (ref_ndof);
			high_rhs->multiply_with_vector(coeff_vec, high_double); 
			vec_high_rhs->zero(); vec_high_rhs->add_vector(high_double);

			UMFPackLinearSolver<double> * highOrd = new UMFPackLinearSolver<double> (high_matrix,vec_high_rhs);	
			if(highOrd->solve()){ 
				u_H = highOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_H, ref_space, &high_sln);	
			  }else error ("Matrix solver failed.\n");

   		/*	sprintf(title, "high-sln adap-step %i", ps);
			  Lowview.set_title(title);
			 Lowview.show(&high_sln);*/ 

		//---------------------------------------antidiffusive fluxes-----------------------------------	

			//	info("assemble fluxes");	
		//	smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
			 antidiffusiveFlux(dof_list,dof_list_2,mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, u_L,coeff_vec, flux_double, 
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);
		
			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);
			vec_rhs->add_vector(flux_double);
			UMFPackLinearSolver<double> * newSol = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(newSol->solve()){ 
				u_new_double = newSol->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_new_double, ref_space, &u_new);	
			}else error ("Matrix solver failed.\n");	 

			if(ps==1) 
					p1_sln.copy(&u_new);		

			 // Visualize the solution.
 		//	sprintf(title, "p1-sln adap-step %i", ps);
			//  Lowview.set_title(title);
			// Lowview.show(&p1_sln);	 
			 
			sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i,ps=%i,", current_time,ts,ps);
				 sview.set_title(title);
					sview.show(&u_new);
				
				//mview.show(ref_space);
	//View::wait(HERMES_WAIT_KEYPRESS);

				int* elements_to_refine = new int[ref_space->get_mesh()->get_max_element_id()+1];   // 0 = nothing..
				int* no_of_refinement_steps = new int[ref_space->get_mesh()->get_max_element_id()+1];	
	
				for(int i = 0; i <= ref_space->get_mesh()->get_max_element_id(); i++){ 
													elements_to_refine[i]= 0;										
															no_of_refinement_steps[i]=0;
					}	
if(ps==2) u_prev_time.copy(&u_new);

smoothness_indicator(ref_space,&u_new,&R_h_1,&R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
			changed = h_p_adap(ref_space,&u_new,&p1_sln,&e_h, &e_h_p,&R_h_1,&R_h_2, &residual, &massmatrix,
			adapting,dof_list,dof_list_2, al, &p1_elements,&neighbor,elements_to_refine,no_of_refinement_steps,h_min,h_max, ts,ps);			
			sprintf(title, "nach changed Mesh, ps=%i, ts=%i", ps,ts);
			mview.set_title(title);
				mview.show(ref_space);

			

			ps++;	
	  // Clean up.
			delete[] flux_double; 
			delete vec_rhs;
			delete[] lumped_double;
			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
			delete lowOrd;
			delete newSol; 	

			delete dof_list;
			delete dof_list_2;		

			delete vec_high_rhs;
			delete highOrd;
			delete [] high_double;


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
			delete [] R_plus;
			delete [] R_minus;
			delete lumped_matrix; 
			delete diffusion;
			delete[] coeff_vec_2;
			if(coeff_vec!=NULL){ delete [] coeff_vec; 	coeff_vec = NULL;}
		  low_matrix->free();
	 		lowmat_rhs->free();

	}while(ps<3);
		

			 // Visualize the solution.
//sprintf(title, "End: Time %3.2f, timestep=%i", current_time,ts);
			 // sview.set_title(title);
			// Lowview.show(&low_sln);	 
		//sview.show(&u_new);
		//mview.show(ref_space);
	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
    // Copy last reference solution into u_prev_time.
  //  u_prev_time.copy(&u_new);



}
while (current_time < T_FINAL);

//lin.save_solution_vtk(&u_prev_time, "end_hpadap_neu.vtk", "solution", mode_3D);
/*sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&u_new);*/

		delete dp_convection;
		delete dp_mass; 
		delete adapting;
		delete ref_space; 

	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;
	delete high_matrix;
	delete high_rhs;

delete al;




  // Wait for the view to be closed.
  View::wait();
  return 0;
}

