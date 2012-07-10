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


const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX = 1; 
const double h_max = 0.1;                       
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.


const double EPS = 1e-8;
const double EPS_smooth = 1e-8;



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
  // Time measurement.
 // TimePeriod cpu_time;
 // cpu_time.tick();

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

 // int ndof = space.get_num_dofs();
 // info("ndof = %d", ndof);

 // Initialize solution of lower & higher order
  Solution<double>  low_sln, u_new,high_sln;
		Solution<double> R_h_1, R_h_2;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);   	



  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, &u_prev_time);
	CustomWeakFormConvection  convection(&u_prev_time);


  // Initialize views.
	//ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(&u_prev_time, HERMES_EPS_HIGH); 

	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
	//mview.show(&space);


  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
//lin.save_solution_vtk(&u_prev_time, "init_h_adap_neu.vtk", "u", mode_3D);

		  // Initialize
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 

	UMFPackMatrix<double> * mass_matrix_uni = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix_uni = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double> * lowmat_rhs_uni = new UMFPackMatrix<double> ; 


		UMFPackMatrix<double> * high_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * high_rhs = new UMFPackMatrix<double> ; 


	double* u_H =NULL;


// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;  //p-adapt-schritte

	 double diag;double f;
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8; 
	info("h_min=%f", h_min);

	AsmList<double>*  al = new AsmList<double>;	

	H1Space<double>* ref_space = new H1Space<double>(&mesh, &bcs, P_INIT);	
	DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, ref_space);
	DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, ref_space);
	HPAdapt * adapting = new HPAdapt(ref_space, HERMES_L2_NORM);

	double* Ax_mass;
	int* Ai_mass;
	int* Ap_mass;
	int ref_ndof;



				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix_uni); 	
			UMFPackMatrix<double> * lumped_matrix_uni = massLumping(mass_matrix_uni);
			//------------------------artificial DIFFUSION D---------------------------------------	
			dp_convection->assemble(conv_matrix_uni, NULL,true);
			UMFPackMatrix<double> * diffusion_uni = artificialDiffusion(conv_matrix_uni);

			lowmat_rhs_uni->create(conv_matrix_uni->get_size(),conv_matrix_uni->get_nnz(), conv_matrix_uni->get_Ap(), conv_matrix_uni->get_Ai(),conv_matrix_uni->get_Ax());
			lowmat_rhs_uni->add_matrix(diffusion_uni); 
			//(1-theta)(K+D)
			if(theta ==1) lowmat_rhs_uni->zero();
			else lowmat_rhs_uni->multiply_with_Scalar((1.0-theta));
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs_uni->add_matrix(lumped_matrix_uni);
			lumped_matrix_uni->multiply_with_Scalar(time_step);  // M_L
			mass_matrix_uni->multiply_with_Scalar(time_step);  // massmatrix = M_C

	int* elements_to_refine = new int[ref_space->get_mesh()->get_max_element_id()];   // 0 = nothing..
	int* no_of_refinement_steps = new int[ref_space->get_mesh()->get_max_element_id()];	
		double* elem_error = new double[ref_space->get_mesh()->get_max_element_id()];


//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 
	
mesh.copy(&basemesh);
	ref_space->set_mesh(&mesh);	ref_space->set_uniform_order(1); ref_space->assign_dofs(); 
	ps=1;

//Adaptivity loop
	do
	{			  // Time measurement.
  //cpu_time.tick(HERMES_SKIP);


			 ref_ndof = ref_space->get_num_dofs();
			dp_mass->set_spaces(ref_space);
			dp_convection->set_spaces(ref_space);

			info(" adap- step %d, timestep %d,ndof = %d ", ps, ts, ref_ndof); 

			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];

			double* lumped_double = new double[ref_ndof];


			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];	
			double* Q_plus_old = new double[ref_ndof]; double* Q_minus_old = new double[ref_ndof];	
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];	

		int* smooth_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[ref_ndof]; 


if(ps==1){

			for(int i=0; i<ref_ndof;i++){Q_plus_old[i]=0.;Q_minus_old[i]=0.;P_plus[i]=0.;}
			//p1_list_fast(ref_space, al, P_plus,&u_prev_time);


			//Initialisierung von Q_plus_old,Q_minus_old
	/*Ax_mass = mass_matrix_uni->get_Ax();
	Ai_mass = mass_matrix_uni->get_Ai();
	Ap_mass = mass_matrix_uni->get_Ap();

		for(int j = 0; j<ref_ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if((Ax_mass[indx]!=0.)&&(j<i)){
						f = lumped_matrix_uni->get_Ax()[i]*(P_plus[j]- P_plus[i])/time_step; 
						if(f>Q_plus_old[i]) Q_plus_old[i] = f;				
						if(f<Q_minus_old[i]) Q_minus_old[i] = f;			
						f= lumped_matrix_uni->get_Ax()[j]*(P_plus[i]- P_plus[j])/time_step; 
						if(f>Q_plus_old[j]) Q_plus_old[j] = f;	
						if(f<Q_minus_old[j]) Q_minus_old[j] = f;
					}
				}
			}*/


	
//----------------- Project the initial condition on the FE space->coeff_vec	
			//info("projection");
			Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix_uni);
			OGProjection<double>::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			//	Solution<double> ::vector_to_solution(coeff_vec, ref_space, &low_sln);
	//	smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2,smooth_elem,smooth_dof,al);
			lumped_flux_limiter(mass_matrix_uni, lumped_matrix_uni, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old,  R_plus, R_minus);

			//Solution<double>::vector_to_solution(coeff_vec, ref_space, &u_new);
		//Solution<double> ::vector_to_solution(coeff_vec_2, ref_space, &high_sln);
	/*	sprintf(title, "proj. FCT_Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(&u_new);
	/*sprintf(title, "proj. lumped_Loesung, ps=%i, ts=%i", ps,ts);
			Lowview.set_title(title);
			Lowview.show(&low_sln);
	sprintf(title, "proj. high_Loesung, ps=%i, ts=%i", ps,ts);
			hview.set_title(title);
			hview.show(&high_sln);*/

//lin.save_solution_vtk(&u_new, "proj_FCT_2.vtk", "u", mode_3D);
	//lin.save_solution_vtk(&low_sln, "lumped_proj.vtk", "u", mode_3D);
		//lin.save_solution_vtk(&high_sln, "high_proj.vtk", "u", mode_3D);
	//View::wait(HERMES_WAIT_KEYPRESS);


	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	
//coeff_vec = u_old = u_n	
			lowmat_rhs_uni->multiply_with_vector(coeff_vec, lumped_double); 

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
				for(int i=0; i<ref_ndof;i++) coeff_vec_2[i]=lumped_double[i]*time_step/lumped_matrix_uni->get_Ax()[i];	
				//coeff_vec_2 = u_L
					Solution<double> ::vector_to_solution(coeff_vec_2, ref_space, &low_sln);	

 			changed = h_p_adap(ref_space,mass_matrix_uni,&low_sln,&R_h_1,&R_h_2,adapting,h_min,h_max,elements_to_refine,no_of_refinement_steps, elem_error);		
			/*sprintf(title, "nach changed Mesh, ps=%i, ts=%i", ps,ts);
			mview.set_title(title);
				mview.show(ref_space);*/



}else{
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ref_ndof);
			double* coeff_vec_3 = new double[ref_ndof]; 
			for(int i=0; i<ref_ndof;i++){Q_plus_old[i]=0.;Q_minus_old[i]=0.;P_plus[i]=0.;}
			//p1_list_fast(ref_space, al, P_plus,&u_prev_time);

				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);
			//------------------------artificial DIFFUSION D---------------------------------------	
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
			
			high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_matrix->multiply_with_Scalar(-theta);
			high_matrix->add_matrix(mass_matrix);  
			high_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_matrix(mass_matrix); 


			//Initialisierung von Q_plus_old,Q_minus_old
	/*Ax_mass = mass_matrix->get_Ax();
	Ai_mass = mass_matrix->get_Ai();
	Ap_mass = mass_matrix->get_Ap();

		for(int j = 0; j<ref_ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if((Ax_mass[indx]!=0.)&&(j<i)){
						f = lumped_matrix->get_Ax()[i]*(P_plus[j]- P_plus[i]); 
						if(f>Q_plus_old[i]) Q_plus_old[i] = f;				
						if(f<Q_minus_old[i]) Q_minus_old[i] = f;			
						f= lumped_matrix->get_Ax()[j]*(P_plus[i]- P_plus[j]); 
						if(f>Q_plus_old[j]) Q_plus_old[j] = f;	
						if(f<Q_minus_old[j]) Q_minus_old[j] = f;
					}
				}
			}*/

			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C

		
//----------------- Project the initial condition on the FE space->coeff_vec	
			//info("projection");
			Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
			OGProjection<double>::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
				Solution<double> ::vector_to_solution(coeff_vec, ref_space, &low_sln);
		smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2,smooth_elem,smooth_dof,al);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old,  R_plus, R_minus,smooth_dof);

			//Solution<double>::vector_to_solution(coeff_vec, ref_space, &u_new);
		//Solution<double> ::vector_to_solution(coeff_vec_2, ref_space, &high_sln);
	/*	sprintf(title, "proj. FCT_Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(&u_new);
	/*sprintf(title, "proj. lumped_Loesung, ps=%i, ts=%i", ps,ts);
			Lowview.set_title(title);
			Lowview.show(&low_sln);
	sprintf(title, "proj. high_Loesung, ps=%i, ts=%i", ps,ts);
			hview.set_title(title);
			hview.show(&high_sln);*/

//lin.save_solution_vtk(&u_new, "proj_FCT_2.vtk", "u", mode_3D);
	//lin.save_solution_vtk(&low_sln, "lumped_proj.vtk", "u", mode_3D);
		//lin.save_solution_vtk(&high_sln, "high_proj.vtk", "u", mode_3D);
	//View::wait(HERMES_WAIT_KEYPRESS);


	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	
//coeff_vec = u_old = u_n	
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
				for(int i=0; i<ref_ndof;i++) coeff_vec_2[i]=lumped_double[i]*time_step/lumped_matrix->get_Ax()[i];	
				//coeff_vec_2 = u_L
					Solution<double> ::vector_to_solution(coeff_vec_2, ref_space, &low_sln);	

		//-------------solution of higher order
			high_rhs->multiply_with_vector(coeff_vec, coeff_vec_3);
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_3);

			UMFPackLinearSolver<double> * highOrd = new UMFPackLinearSolver<double> (high_matrix,vec_rhs);	
			if(highOrd->solve()){ 
				u_H = highOrd->get_sln_vector();  
				//Solution<double> ::vector_to_solution(u_H, ref_space, &high_sln);	
			  }else error ("Matrix solver failed.\n");

   		/*	sprintf(title, "high-sln adap-step %i", ps);
			  Lowview.set_title(title);
			 Lowview.show(&high_sln);*/ 

  // CPU time 
  //double time1 = cpu_time.tick().last();
  // Time measurement.
 // cpu_time.tick(HERMES_SKIP);

		//---------------------------------------antidiffusive fluxes-----------------------------------	
			//	info("assemble fluxes");	
			smoothness_indicator(ref_space,&low_sln,&R_h_1,&R_h_2, smooth_elem,smooth_dof,al);
			 antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, coeff_vec_2,coeff_vec, coeff_vec_3, 
									P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old,  R_plus, R_minus,smooth_dof);
		//coeff_vec = u_old, coeff_vec_2 = u_L,coeff_vec_3= flux

			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);
			vec_rhs->add_vector(coeff_vec_3);
			UMFPackLinearSolver<double> * newSol = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(newSol->solve()){ 
				Solution<double> ::vector_to_solution(newSol->get_sln_vector(), ref_space, &u_new);	
			}else error ("Matrix solver failed.\n");	 

  // Measure the projection time.
  //double flux_time = cpu_time.tick().last();
//info("CPU_Time:  %g -------------- %g", time1, flux_time);

			 // Visualize the solution.			 
		/*sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i,ps=%i,", current_time,ts,ps);
				 sview.set_title(title);
					sview.show(&u_new);*/
				//mview.show(ref_space);
	//View::wait(HERMES_WAIT_KEYPRESS);
			
 			u_prev_time.copy(&u_new); 
			delete newSol; 
			delete highOrd;
			high_rhs->free();
			high_matrix->free();
			delete lumped_matrix; 
			delete diffusion;
			delete vec_rhs;
		  low_matrix->free();
	 		lowmat_rhs->free();
			delete[] coeff_vec_3;


		}			

			ps++;	
	  // Clean up.

			delete[] lumped_double;
		delete [] smooth_elem;
		delete [] smooth_dof;
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] Q_plus_old;
			delete [] Q_minus_old;
			delete [] R_plus;
			delete [] R_minus;
			delete[] coeff_vec_2;
			 delete [] coeff_vec;



	}while(ps<3);
		


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

}
while (current_time < T_FINAL);

CustomInitialCondition exact_solution(ref_space->get_mesh());
	Adapt<double>* error_estimation = new Adapt<double>(ref_space, HERMES_L2_NORM);	
double err_est = error_estimation->calc_err_est(&exact_solution,&u_prev_time,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
double err_est_2 = error_estimation->calc_err_est(&u_prev_time,&exact_solution,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
printf("err_est = %f, err_est_2 =%f,  ndof = %d", err_est,err_est_2, ref_ndof);
FILE * pFile;
pFile = fopen ("error.txt","w");
     fprintf (pFile, "err_est = %f, err_est_2 =%f,  ndof = %d", err_est,err_est_2, ref_ndof);
fclose (pFile);

lin.save_solution_vtk(&u_prev_time, "end_non_smooth_h_adap.vtk", "solution", mode_3D);


//sview.show(&u_new);
//mview.show(ref_space);
//mview.save_numbered_screenshot("solution.bmp", true);

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

	delete mass_matrix_uni;  
	delete conv_matrix_uni;
	delete lowmat_rhs_uni;
delete al;


			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
		delete [] elem_error;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

