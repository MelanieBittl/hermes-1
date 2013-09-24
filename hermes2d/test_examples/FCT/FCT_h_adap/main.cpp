#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "prev_solution.h"
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
const int P_MAX = 3; 
const double h_max = 0.1;                       
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.
 
const double EPS = 1e-10;

const double THRESHOLD = 0.3;

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
#include "reg_estimator.cpp"
#include "z_z.cpp"
#include "error_estimations.cpp"


int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);


  // Initialize boundary conditions.
  DefaultEssentialBCConst<double>  bc_essential(BDY_IN, 0.0);
  EssentialBCs<double>  bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
//SpaceSharedPtr<double> space(new SpaceBB<double>(mesh, P_INIT));	

	SpaceSharedPtr<double> space(new H1Space<double>(mesh,&bcs, P_INIT));	

 // Initialize solution of lower & higher order
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double>  low_sln(new Solution<double>);
  MeshFunctionSharedPtr<double>  high_sln(new Solution<double>);
  MeshFunctionSharedPtr<double>  u_prev(new PrevSolution);
  MeshFunctionSharedPtr<double>  R_h_1(new Solution<double>);
  MeshFunctionSharedPtr<double>  R_h_2(new Solution<double>);


  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));	
	
  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, u_prev_time);
	CustomWeakFormConvection  convection(u_prev_time);

  // Initialize views.
	//ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(u_prev_time);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
	//mview.show(space);


  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
//lin.save_solution_vtk(u_prev_time, "init_h_adap_neu.vtk", "u", mode_3D);

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


	UMFPackMatrix<double> * mass_matrix_uni = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix_uni = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double> * lowmat_rhs_uni = new UMFPackMatrix<double> ; 
			DiscreteProblem<double> * dp_mass_uni = new DiscreteProblem<double> (&massmatrix, space);
			DiscreteProblem<double> * dp_convection_uni = new DiscreteProblem<double> (&convection, space);
				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass_uni->assemble(mass_matrix_uni); 	
			UMFPackMatrix<double> * lumped_matrix_uni = massLumping(mass_matrix_uni);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection_uni->assemble(conv_matrix_uni);
			UMFPackMatrix<double> * diffusion_uni = artificialDiffusion(conv_matrix_uni);
			//--------------------------------------------------------------------------------------------
			lowmat_rhs_uni->create(conv_matrix_uni->get_size(),conv_matrix_uni->get_nnz(), conv_matrix_uni->get_Ap(), conv_matrix_uni->get_Ai(),conv_matrix_uni->get_Ax());
			lowmat_rhs_uni->add_matrix(diffusion_uni); 
			if(theta ==1) lowmat_rhs_uni->zero();
			else lowmat_rhs_uni->multiply_with_Scalar((1.0-theta)); 
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs_uni->add_matrix(lumped_matrix_uni);			
			lumped_matrix_uni->multiply_with_Scalar(time_step);  // M_L
			mass_matrix_uni->multiply_with_Scalar(time_step);  // massmatrix = M_C


// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;  //p-adapt-schritte

	int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   // 0 = nothing..
	int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];	
		double* elem_error = new double[space->get_mesh()->get_max_element_id()];

	 double diag,f;
	Element* e =NULL;
	for_all_active_elements(e, mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8; 


	AsmList<double> al;

			DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, space);
			DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, space);

		int ndof;
    OGProjection<double> ogProjection;
  DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
  AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);

//Timestep loop
do
{	  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time); 	 
	
 	 mesh->copy(basemesh);
	space->set_mesh(mesh);	space->set_uniform_order(P_INIT); space->assign_dofs(); 
	ps=1; 
dp_mass->set_space(space); 
dp_convection->set_space(space); 
//Adaptivity loop
	do
	{		
			ndof = space->get_num_dofs();
  	Hermes::Mixins::Loggable::Static::info(" adap- step %i, ndof = %i  ", ps, ndof); 	 

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];

			double* lumped_double = new double[ndof];
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	
		int* smooth_elem = new int[space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[ndof];
 
			for(int i=0; i<ndof;i++){Q_plus_old[i]=0.;Q_minus_old[i]=0.;P_plus[i]=0.;}

if(ps==1){
			HPAdapt * adapting = new HPAdapt(space,&error_calculator);
	//--------Project the initial condition on the FE space->coeff_vec	---------------
			if(ts==1) {
				Lumped_Projection::project_lumped(space, u_prev_time, coeff_vec,   lumped_matrix_uni);
    ogProjection.project_global(space,u_prev_time, coeff_vec_2,  HERMES_L2_NORM);
			}else{
				 Lumped_Projection::project_lumped(space, u_prev, coeff_vec,  lumped_matrix_uni);
    ogProjection.project_global(space,u_prev, coeff_vec_2,  HERMES_L2_NORM);

			}
			Solution<double>::vector_to_solution(coeff_vec, space, low_sln);
//		smoothness_indicator(space,low_sln,R_h_1,R_h_2,smooth_elem,smooth_dof,&al,mass_matrix);

			lumped_flux_limiter(mass_matrix_uni, lumped_matrix_uni, coeff_vec, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus);

			//Solution<double>::vector_to_solution(coeff_vec, space, u_new);

			/*sprintf(title, "proj. Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(u_new);*/

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	
//coeff_vec = u_old = u_n	
			lowmat_rhs_uni->multiply_with_vector(coeff_vec, lumped_double); 

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
				for(int i=0; i<ndof;i++) coeff_vec_2[i]=lumped_double[i]*time_step/lumped_matrix_uni->get_Ax()[i];	
				//coeff_vec_2 = u_L
					Solution<double> ::vector_to_solution(coeff_vec_2, space, low_sln);
			changed = h_p_adap(space,mass_matrix_uni,low_sln,R_h_1,R_h_2,adapting,h_min,h_max,elements_to_refine, no_of_refinement_steps,elem_error);	

			//sprintf(title, "nach changed Mesh, ps=%i, ts=%i", ps,ts);
			//mview.set_title(title);
			//mview.show(space);		
			//View::wait(HERMES_WAIT_KEYPRESS);
		delete adapting;

}else{
			double* coeff_vec_3 = new double[ndof]; 

				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection->assemble(conv_matrix);
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
			UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ndof);
		high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_matrix->multiply_with_Scalar(-theta);
			high_matrix->add_matrix(mass_matrix);  
			high_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_matrix(mass_matrix); 	

			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C

	//--------Project the initial condition on the FE space->coeff_vec	---------------
			if(ts==1) {
				Lumped_Projection::project_lumped(space, u_prev_time, coeff_vec,lumped_matrix);
    		ogProjection.project_global(space,u_prev_time, coeff_vec_2,  HERMES_L2_NORM);
			}else{
				 Lumped_Projection::project_lumped(space, u_prev, coeff_vec,lumped_matrix);
		    ogProjection.project_global(space,u_prev, coeff_vec_2,  HERMES_L2_NORM);
			}
		/*	Solution<double>::vector_to_solution(coeff_vec, space, low_sln);
		smoothness_indicator(space,low_sln,R_h_1,R_h_2,smooth_elem,smooth_dof,&al,mass_matrix);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus,smooth_dof);*/

lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus);

			/*Solution<double>::vector_to_solution(coeff_vec, space, u_new);
			sprintf(title, "proj. Loesung, ps=%i, ts=%i", ps,ts);
			pview.set_title(title);
			pview.show(u_new);

			mview.show(space);		
			View::wait(HERMES_WAIT_KEYPRESS);*/

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	
//coeff_vec = u_old = u_n	
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
				for(int i=0; i<ndof;i++) coeff_vec_2[i]=lumped_double[i]*time_step/lumped_matrix->get_Ax()[i];	
				//coeff_vec_2 = u_L

		//-------------solution of higher order
			high_rhs->multiply_with_vector(coeff_vec, coeff_vec_3); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_3);
			UMFPackLinearMatrixSolver<double> * highOrd = new UMFPackLinearMatrixSolver<double> (high_matrix,vec_rhs);	
			  try
  {
   highOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_H = highOrd->get_sln_vector();  
				//Solution<double> ::vector_to_solution(u_H, space, high_sln);	


		//---------------------------------------antidiffusive fluxes-----------------------------------	
	//	smoothness_indicator(space,low_sln,R_h_1,R_h_2,smooth_elem,smooth_dof,&al,mass_matrix);
//antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, coeff_vec_2,coeff_vec, coeff_vec_3, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus,smooth_dof);
antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, coeff_vec_2,coeff_vec, coeff_vec_3, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus);

			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);
			vec_rhs->add_vector(coeff_vec_3);
			UMFPackLinearMatrixSolver<double> * newSol = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
  try
  {
   newSol->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_new_double = newSol->get_sln_vector();  
				Solution<double>::vector_to_solution(u_new_double, space, u_new);	


			 // Visualize the solution.		 
				sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i,ps=%i,", current_time,ts,ps);
				 sview.set_title(title);
mview.show(space);
					sview.show(u_new);

		u_prev->copy(u_new);
		u_prev->set_own_mesh(space->get_mesh());


		delete newSol;
		delete highOrd;
			delete vec_rhs;
			high_matrix->free();
			high_rhs->free();
		  low_matrix->free();
		delete[] coeff_vec_3; 
			delete lumped_matrix; 
			delete diffusion;

}//End ps=2			


			ps++;	
	  // Clean up.
			delete[] lumped_double;
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
	 		lowmat_rhs->free();
		delete [] smooth_elem;
		delete [] smooth_dof;





	}while(ps<3);		

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
 	ts++;

}
while (current_time < T_FINAL);

//lin.save_solution_vtk(u_prev, "end_h_adap_neu.vtk", "solution", mode_3D);
/*sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(u_new);*/

calc_error_total(u_prev, u_prev_time,space);


		delete dp_convection;
		delete dp_mass; 

	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;
	delete high_matrix;
	delete high_rhs;

	delete mass_matrix_uni;  
	delete conv_matrix_uni;
	delete lowmat_rhs_uni;
		delete dp_convection_uni;
		delete dp_mass_uni; 
			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
		delete [] elem_error;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

