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
using namespace Hermes::Solvers;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 


const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX = 2; 
const double h_max = 0.1;                       
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.
 
const double THRESHOLD = 0.3;
 
const double EPS = 1e-10;
const double EPS_smooth = 1e-10;



const int NDOF_STOP = 20000;   

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

//FCT & p-Adaptivity
#include "fct.cpp"
#include "p_adapt.cpp"
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "p1_list.cpp"
#include "reg_estimator.cpp"
#include "z_z.cpp"



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
SpaceSharedPtr<double> space(new H1Space<double>(mesh,&bcs, P_INIT));	

  int ndof = space->get_num_dofs();
  

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
	//Lowview.show(u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
//ScalarView hview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));



  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
Orderizer ord;


		  // Initialize
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;   //M_c/tau
	CSCMatrix<double> * conv_matrix = new CSCMatrix<double> ;   //K
	CSCMatrix<double> * low_matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * lowmat_rhs = new CSCMatrix<double> ; 


		CSCMatrix<double> * high_matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * high_rhs = new CSCMatrix<double> ; 

	double* u_L = NULL; 
	double* u_new_double = NULL; 
	double* u_H =NULL;
		double* u_proj_fct =NULL;




// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int as = 1;  //p-adapt-schritte


	 double diag;			double f;
	Element* e =NULL;
	for_all_active_elements(e, mesh){diag = e->get_diameter(); break;}
	const double h_start = diag;
	const double h_min = diag/8; 


	AsmList<double>*  al = new AsmList<double>;	

  OGProjection<double> ogProjection;
  DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
  AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);

	SpaceSharedPtr<double> ref_space(new H1Space<double>(mesh,&bcs, P_INIT));	

//Timestep loop
do
{	 

		Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time);
  
	  
	if(ts>1){
			mesh->copy(basemesh);
			ref_space->set_mesh(mesh);
			ref_space->set_uniform_order(P_INIT); ref_space->assign_dofs(); //mview.show(ref_space);
	}
	

	as=1; 


	do
	{	
			int ref_ndof = ref_space->get_num_dofs(); 
			Hermes::Mixins::Loggable::Static::info(" adap- step %d, timestep %d,ndof = %d ", as, ts, ref_ndof);

		DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, ref_space);
		DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, ref_space);
		HPAdapt * adapting = new HPAdapt(ref_space,&error_calculator);

			double* coeff_vec = new double[ref_ndof];
			double* coeff_vec_2 = new double[ref_ndof];
			double* coeff_vec_3 = new double[ref_ndof];
			double* P_plus = new double[ref_ndof]; double* P_minus = new double[ref_ndof];
			double* Q_plus = new double[ref_ndof]; double* Q_minus = new double[ref_ndof];		
			double* R_plus = new double[ref_ndof]; double* R_minus = new double[ref_ndof];

		int* smooth_elem = new int[ref_space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[ref_ndof];
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ref_ndof);
			double* lumped_double = new double[ref_ndof];



		if(as==1){

		//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix,vec_rhs); 	
			CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection->assemble(conv_matrix);
			CSCMatrix<double> * diffusion = artificialDiffusion(conv_matrix);
			//--------------------------------------------------------------------------------------------

			lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			lowmat_rhs->add_sparse_matrix(diffusion); 
			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			//(-theta)(K+D)
			if(theta==0) low_matrix->zero();
			else	low_matrix->multiply_with_Scalar(-theta);
			//(1-theta)(K+D)
			if(theta ==1) lowmat_rhs->zero();
			else lowmat_rhs->multiply_with_Scalar((1.0-theta));

			//M_L/tau - theta(D+K)
			low_matrix->add_sparse_matrix(lumped_matrix);  
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs->add_sparse_matrix(lumped_matrix);	



			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C

			//------------------ Project the initial condition on the FE space->coeff_vec	--------------
			if(ts==1) {
				Lumped_Projection::project_lumped(ref_space, u_prev_time, coeff_vec, lumped_matrix);
    ogProjection.project_global(ref_space,u_prev_time, coeff_vec_2,  HERMES_L2_NORM);
			}else{
				 Lumped_Projection::project_lumped(ref_space, u_prev, coeff_vec,lumped_matrix);
    ogProjection.project_global(ref_space,u_prev, coeff_vec_2,  HERMES_L2_NORM);

			}
			Solution<double>::vector_to_solution(coeff_vec, ref_space, low_sln);
		smoothness_indicator(ref_space,low_sln,R_h_1,R_h_2,smooth_elem,smooth_dof,al,mass_matrix,true);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus,smooth_dof);



			Solution<double>::vector_to_solution(coeff_vec, ref_space, u_new);
		sprintf(title, "proj. Loesung, as=%i, ts=%i", as,ts);
			pview.set_title(title);
			pview.show(u_new);
//View::wait(HERMES_WAIT_KEYPRESS);

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 

	//-------------------------solution of lower order  M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n------------	
				for(int i=0; i<ref_ndof;i++) coeff_vec_2[i]=lumped_double[i]*time_step/lumped_matrix->get_Ax()[i];	
				// u_L = coeff_vec_2 
					Solution<double> ::vector_to_solution(coeff_vec_2, ref_space, low_sln);	

				//smoothness_indicator(ref_space,low_sln,R_h_1,R_h_2, smooth_elem,smooth_dof,al,mass_matrix,true);
			smoothness_indicator(ref_space,low_sln,R_h_1,R_h_2, smooth_elem,smooth_dof,al,mass_matrix,true);
			changed = h_p_adap(ref_space,u_prev_time, low_sln,R_h_1,R_h_2,&massmatrix, adapting,al, h_min,h_max, ts,as,smooth_elem, h_start);	
			
		
sprintf(title, "nach changed Mesh, as=%i, ts=%i", as,ts);
			mview.set_title(title);
				mview.show(ref_space);

View::wait(HERMES_WAIT_KEYPRESS);
			delete lumped_matrix; 
			delete diffusion;
	


			}else{


//P=2 -----------------------------( nicht uniformes Gitter!!!!!!)
					bool* fct = new bool[ref_ndof]; 
			for(int i=0; i<ref_ndof;i++){fct[i]=false; /*P_plus[i]=0.;*/ }
					//if(ts==1) p1_list(ref_space, fct, al,P_plus,u_prev_time,h_start);		
					//else p1_list(ref_space, fct, al,P_plus,u_prev,h_start);		
					p1_list(ref_space, fct, al,h_start);		


				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix,vec_rhs); 	
				CSCMatrix<double>* lumped_matrix = massLumping(fct,mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection->assemble(conv_matrix);
				CSCMatrix<double>* diffusion = artificialDiffusion(fct,conv_matrix);

			//--------------------------------------------------------------------------------------------
			lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			lowmat_rhs->add_sparse_matrix(diffusion); 
			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			//(-theta)(K+D)
			if(theta==0) low_matrix->zero();
			else	low_matrix->multiply_with_Scalar(-theta);
			//(1-theta)(K+D)
			if(theta ==1) lowmat_rhs->zero();
			else lowmat_rhs->multiply_with_Scalar((1.0-theta));

			//M_L/tau - theta(D+K)
			low_matrix->add_sparse_matrix(lumped_matrix);  
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs->add_sparse_matrix(lumped_matrix);	


			high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_matrix->multiply_with_Scalar(-theta);
			high_matrix->add_sparse_matrix(mass_matrix);  
			high_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_sparse_matrix(mass_matrix); 


			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C

			// Project the initial condition on the FE space->coeff_vec	
			if(ts==1) {
				Lumped_Projection::project_lumped(ref_space, u_prev_time, coeff_vec,lumped_matrix);
   		  ogProjection.project_global(ref_space,u_prev_time, coeff_vec_2,  HERMES_L2_NORM);

			}else{
				 Lumped_Projection::project_lumped(ref_space, u_prev, coeff_vec, lumped_matrix);
    		 ogProjection.project_global(ref_space,u_prev, coeff_vec_2,  HERMES_L2_NORM);
			}

			Solution<double>::vector_to_solution(coeff_vec, ref_space, low_sln);
			smoothness_indicator(ref_space,low_sln,R_h_1,R_h_2,smooth_elem,smooth_dof,al,mass_matrix,true);
				lumped_flux_limiter(fct,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
										P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus,smooth_dof);


		/*	Solution<double>::vector_to_solution(coeff_vec, ref_space, u_new);
			sprintf(title, "proj. Loesung, as=%i, ts=%i", as,ts);
			pview.set_title(title);
			pview.show(u_new);*/

	
		//	Hermes::Mixins::Loggable::Static::info(" Solve linear system ");

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 
			vec_rhs->zero(); vec_rhs->add_vector(lumped_double);

	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
			lumped_matrix->multiply_with_Scalar(1./time_step); //M_L/tau
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (lumped_matrix,vec_rhs);	
  try
  {
   lowOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_L = lowOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_L, ref_space, low_sln);	

			lumped_matrix->multiply_with_Scalar(time_step);  // M_L


		//-------------solution of higher order------
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
				Solution<double> ::vector_to_solution(u_H, ref_space, high_sln);	

   		/*	sprintf(title, "high-sln adap-step %i", as);
			  hview.set_title(title);
			 hview.show(high_sln);*/

		//---------------------------------------antidiffusive fluxes-----------------------------------		
			smoothness_indicator(ref_space,low_sln,R_h_1,R_h_2, smooth_elem,smooth_dof,al,mass_matrix,true);
			 antidiffusiveFlux(fct,mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, u_L,coeff_vec, coeff_vec_3, 
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus,smooth_dof);
	
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
				Solution<double> ::vector_to_solution(newSol->get_sln_vector(), ref_space, u_new);	



			 // Visualize the solution.		 
			sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
				 sview.set_title(title);
					sview.show(u_new);
				
				//mview.show(ref_space);
	//View::wait(HERMES_WAIT_KEYPRESS);



			 u_prev->copy(u_new); 
			 u_prev->set_own_mesh(u_new->get_mesh());
			 
			delete lumped_matrix; 
			delete diffusion;
			delete highOrd;
			delete lowOrd;
			delete newSol; 	
			delete [] fct;
			high_rhs->free();
			high_matrix->free();

}//Ende as hoeher
			

			as++;	
	  // Clean up.

		delete[] coeff_vec_3; 
		delete vec_rhs;
		delete[] lumped_double;


		delete [] smooth_elem;
		delete [] smooth_dof;
		delete [] P_plus;
		delete [] P_minus;
		delete [] Q_plus;
		delete [] Q_minus;
		delete [] R_plus;
		delete [] R_minus;
		delete[] coeff_vec_2;
		delete [] coeff_vec; 
		low_matrix->free();
		lowmat_rhs->free();

		delete dp_convection;
		delete dp_mass; 
		delete adapting;
		

	}while(as<3);

			 // Visualize the solution.
//sprintf(title, "End: Time %3.2f, timestep=%i", current_time,ts);
			//  sview.set_title(title);
			// Lowview.show(low_sln);	 
		//sview.show(u_new);
		//mview.show(ref_space);
	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
  


}
while (current_time < T_FINAL);
Hermes::Mixins::Loggable::Static::info(" END");


lin.save_solution_vtk(u_prev, "end_hpadap_smooth.vtk", "solution", mode_3D);
ord.save_mesh_vtk(ref_space, "end_mesh");
ord.save_orders_vtk(ref_space, "end_order.vtk");




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

