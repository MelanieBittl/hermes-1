#define HERMES_REPORT_ALL
#include "definitions.h"



using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 


const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
                    
const double time_step = 0.01;                           // Time step.

const double T_FINAL = 1000;                       // Time interval length.



const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 20;                  // Maximum allowed number of Newton iterations.

const int NDOF_STOP = 20000;   


MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


int main(int argc, char* argv[])
{
   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, P_INIT);

  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

 // Initialize solution of lower & higher order
  Solution<double>  high_sln;

  // Previous time level solution (initialized by the initial condition).
		ConstantSolution<double> u_prev_time(&mesh, 0.0);	

	CustomExactSolution exact_solution(&mesh);

	H1Space<double>* ref_space = new H1Space<double>(&mesh, P_INIT);	

	ConvectionForm high_form(time_step, &u_prev_time, &exact_solution);


  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;

	DiscreteProblem<double> * dp = new DiscreteProblem<double> (&high_form, ref_space);
	double* u_H =NULL;
	int ref_ndof = ref_space->get_num_dofs();
	double* coeff_vec = new double[ref_ndof];
  memset(coeff_vec, 0, ref_ndof*sizeof(double));

		UMFPackMatrix<double> * high_matrix = new UMFPackMatrix<double> ;
		UMFPackVector<double> * high_rhs = new UMFPackVector<double> (ref_ndof); 
//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

		dp->assemble(high_matrix, high_rhs);  //M+K
		UMFPackLinearSolver<double> * highOrd = new UMFPackLinearSolver<double> (high_matrix,high_rhs);	
			if(highOrd->solve()){ 
				u_H = highOrd->get_sln_vector();  
				Solution<double> ::vector_to_solution(u_H, ref_space, &u_prev_time);	
			  }else error ("Matrix solver failed.\n");



   		sprintf(title, "high-sln time-step %i", ts);
			  Lowview.set_title(title);
			 Lowview.show(&u_prev_time);	



	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;


}
while (current_time < T_FINAL);


delete dp;

delete high_matrix;
delete high_rhs;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

