#define HERMES_REPORT_ALL
#include "definitions.h"


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM = 4;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.
                     
const double time_step = 1e-3;                           // Time step.
const double T_FINAL = 2*PI;                         // Time interval length. 
//const double T_FINAL = 0.5*PI; 


const double theta = 1.;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

const bool all = true;
const bool DG = true;
const bool SD = false;

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";


#include "error_estimations.cpp"


int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);
 /*  MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(basemesh);
   View::wait();*/

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);


  // Initialize boundary conditions.
//CustomDirichletCondition bc_essential(Hermes::vector<std::string>("inlet1","inlet2"));
 // EssentialBCs<double>  bcs(&bc_essential);
  
  // Create an space with default shapeset.  

  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));	

  int ndof = space->get_num_dofs();
  
  /* BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
  bview.show(&space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));
 
  // Initialize views.
	ScalarView sview("Loesung", new WinGeom(500, 500, 500, 400));
	ScalarView lview("init Loesung", new WinGeom(500, 0, 500, 400));
	//lview.show(u_prev_time);

  OGProjection<double> ogProjection;
  
	char title[100];
    int ts = 1;
    double current_time =0.;
	
	int ref_ndof = space->get_num_dofs();	
	//Hermes::Mixins::Loggable::Static::info(" ndof = %d ", ref_ndof); 

	Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads,2);
CustomWeakForm wf(u_prev_time, mesh,time_step, theta, all, DG, SD);

 // Initialize linear solver.
 Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  do
  {Hermes::Mixins::Loggable::Static::info("time=%f, ndof = %i ", current_time,ref_ndof); 
					// Solve the linear problem.
					try
					{
						linear_solver.solve();

						// Get the solution vector.
						double* sln_vector = linear_solver.get_sln_vector();

						// Translate the solution vector into the previously initialized Solution.
						Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, u_new);
						//sview.show(u_new);
 						if(ts==1) wf.set_ext(u_new);
					}catch(std::exception& e)
					{
						std::cout << e.what();
					}
	 	current_time += time_step;
	 	ts++;
	  }
  while (current_time<T_FINAL);
  

calc_error_total(u_new, u_prev_time,space);


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

