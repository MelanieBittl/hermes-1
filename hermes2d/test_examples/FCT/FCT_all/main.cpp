#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_NUM =5;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.                      
const double time_step =1e-5;                           // Time step.
const double T_FINAL = 0.5; //2.;                     // Time interval length. 




const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


//FCT & p-Adaptivity
#include "fct.cpp"
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"

#include "error_estimations.cpp"

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

  // Create an H1 space with default shapeset.
SpaceSharedPtr<double> space(new H1Space<double>(mesh,P_INIT));	



 // Initialize solution of lower & higher order
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double>  low_sln(new Solution<double>);
  MeshFunctionSharedPtr<double>  high_sln(new Solution<double>);

  MeshFunctionSharedPtr<double>  u_proj(new Solution<double>);

  MeshFunctionSharedPtr<double>  R_h_1(new Solution<double>);
  MeshFunctionSharedPtr<double>  R_h_2(new Solution<double>);

  // Previous time level solution (initialized by the initial condition).
double current_time = 0.0; 
  MeshFunctionSharedPtr<double> u_init(new CustomInitialCondition(mesh,current_time));	
  MeshFunctionSharedPtr<double> u_exact(new CustomInitialCondition(mesh, current_time));

  // Initialize the weak formulation.
	CustomWeakFormMassmatrix  massmatrix(time_step, u_init);
	CustomWeakFormConvection  convection(true,false,false);
CustomWeakFormConvection  diff_wf(false,true,false);
CustomWeakFormConvection  reac_wf(false,false,true);

CustomWeakFormSurface surface(u_init);
CustomWeakFormRHS wf_rhs(u_exact,u_init,mesh,time_step,theta);

			DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, space);
			DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, space);
			DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&surface, space);
DiscreteProblem<double> * dp_rhs = new DiscreteProblem<double> (&wf_rhs, space);

DiscreteProblem<double> * dp_diff = new DiscreteProblem<double> (&diff_wf, space);
DiscreteProblem<double> * dp_reac = new DiscreteProblem<double> (&reac_wf, space);



  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	MeshView mview("mesh", new WinGeom(0, 0, 500, 400));
OrderView oview("order", new WinGeom(500,500,500,400));
sview.set_min_max_range(-1, 1.);

		  // Initialize
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;   //M_c/tau
	CSCMatrix<double> * conv_matrix = new CSCMatrix<double> ;   //K-L
CSCMatrix<double> * diff_matrix = new CSCMatrix<double> ;
CSCMatrix<double> * reac_matrix = new CSCMatrix<double> ;
	CSCMatrix<double> * low_matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * lowmat_rhs = new CSCMatrix<double> ; 

		CSCMatrix<double> * high_matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * high_rhs = new CSCMatrix<double> ; 

	CSCMatrix<double> * mat_surface = new CSCMatrix<double> ; 

	double* u_H =NULL;	



// Time stepping loop:
	
	int ts = 1;
	char title[100];


		int ndof = space->get_num_dofs();
    OGProjection<double> ogProjection;
int time_counter =1.;


	


			SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
			SimpleVector<double> * vec_surf = new SimpleVector<double> (ndof);
SimpleVector<double> * rhs = new SimpleVector<double> (ndof); 
		double* coeff_vec = new double[ndof]; 
		double* coeff_vec_2 = new double[ndof]; 
		double* coeff_vec_3 = new double[ndof]; 

			double* lumped_double = new double[ndof];

			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	


				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			dp_mass->assemble(mass_matrix); 	
			CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			dp_convection->assemble(conv_matrix);
		CSCMatrix<double> * diffusion = artificialDiffusion(conv_matrix);

dp_diff->assemble(diff_matrix);
conv_matrix->add_sparse_matrix(diff_matrix);

dp_reac->assemble(reac_matrix);
//conv_matrix->add_sparse_matrix(reac_matrix);


CSCMatrix<double> * diffusion_diff = artificialDiffusion(diff_matrix);
CSCMatrix<double> * diffusion_reac = massLumping(reac_matrix);

diffusion->add_sparse_matrix(diffusion_diff);  //diffusion->add_sparse_matrix(diffusion_reac);

			//dp_surf->assemble(mat_surface, vec_surf);
			//--------------------------------------------------------------------------------------------

			lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			lowmat_rhs->add_sparse_matrix(diffusion);
			lowmat_rhs->add_sparse_matrix(diffusion_reac);  

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

	//low_matrix->add_sparse_matrix(mat_surface); 
conv_matrix->add_sparse_matrix(reac_matrix);
			high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_matrix->multiply_with_Scalar(-theta);
			high_matrix->add_sparse_matrix(mass_matrix); 


			high_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_sparse_matrix(mass_matrix); 	

			//high_matrix->add_sparse_matrix(mat_surface); 



		lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C


	//--------Project the initial condition on the FE space->coeff_vec	---------------

	Lumped_Projection::project_lumped(space, u_init, coeff_vec,lumped_matrix);
    		ogProjection.project_global(space,u_init, coeff_vec_2,  HERMES_L2_NORM);
			Solution<double>::vector_to_solution(coeff_vec_2, space, u_proj);

			for(int i=0; i< ndof;i++){Q_plus_old[i]=0.;Q_minus_old[i]=0.;}

lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus);

//Timestep loop
do
{	 
wf_rhs.set_current_time(current_time);
	dynamic_cast<CustomInitialCondition*>(u_exact.get())->set_time(current_time);
	dp_rhs->assemble(rhs);

  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time); 	 
	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	
//coeff_vec = u_old = u_n	
	lowmat_rhs->multiply_with_vector(coeff_vec, lumped_double); 

			for(int i=0; i<ndof;i++)lumped_double[i]+=rhs->get(i);
					//lumped_double[i]+=vec_surf->get(i)+rhs->get(i);
	//-------------------------solution of lower order------------	
			  // Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
				for(int i=0; i<ndof;i++)
			 		coeff_vec_2[i]=lumped_double[i]*time_step/(lumped_matrix->get_Ax()[i]);
				
				//Solution<double> ::vector_to_solution(coeff_vec_2, space, low_sln);	


		//-------------solution of higher order
			high_rhs->multiply_with_vector(coeff_vec, coeff_vec_3); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_3);
//vec_rhs->add_vector(vec_surf);
vec_rhs->add_vector(rhs);
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

antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion,u_H, coeff_vec_2,coeff_vec, coeff_vec_3, P_plus, P_minus, Q_plus, Q_minus,Q_plus_old, Q_minus_old, R_plus, R_minus);

		//mview.show(space->get_mesh());//press m for marker


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

	for(int i =0; i<ndof; i++)	coeff_vec[i] =newSol->get_sln_vector()[i]; 	
				//coeff_vec[i] =coeff_vec_2[i];
 
				Solution<double>::vector_to_solution(coeff_vec, space, u_new);	


			 // Visualize the solution.		
				//sprintf(title, "korrigierte Loesung: Time %3.2f,timestep %i", current_time,ts);
				// sview.set_title(title);
sview.set_min_max_range(0.,0.5);
					sview.show(u_new);

		delete newSol;
		delete highOrd;
	

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
 	ts++;



}
while (current_time < T_FINAL);

calc_error_total(u_new, u_init,space);



	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;
	delete high_matrix;
	delete high_rhs;
delete vec_rhs;
	delete mat_surface; 

delete diff_matrix;
delete reac_matrix;

			delete dp_mass;
			delete dp_convection;
			delete dp_surf;
delete dp_rhs;
delete dp_diff;
delete dp_reac;
delete lumped_matrix;
delete diffusion;

			 delete [] coeff_vec; 
			delete[] coeff_vec_2;
			 delete [] coeff_vec_3; 
			delete[] lumped_double;
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] Q_plus_old;
			delete [] Q_minus_old;
			delete [] R_plus;
			delete [] R_minus;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

