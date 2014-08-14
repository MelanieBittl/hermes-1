#define HERMES_REPORT_ALL
#include "definitions.h"


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.
                     
const double time_step =25e-5;                           // Time step.
const double T_FINAL = 0.5;                         // Time interval length. 


const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)
const double theta_DG =0.5;

const bool all = true;
const bool DG = false;


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
mloader.load("unit.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);
  
  // Create an space with default shapeset.  
 //SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
 SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));	

  int ndof = space->get_num_dofs();

double current_time = 0.;

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_exact(new CustomInitialCondition(mesh, current_time));
 
  // Initialize views.
	ScalarView sview("Loesung", new WinGeom(500, 500, 500, 400));
	ScalarView lview("init Loesung", new WinGeom(500, 0, 500, 400));
	//lview.show(u_exact);
//View::wait(HERMES_WAIT_KEYPRESS);

 
	char title[100];
    int ts = 1;
	
	int ref_ndof = space->get_num_dofs();	

	dynamic_cast<CustomInitialCondition*>(u_exact.get())->set_time(current_time);

	CustomWeakForm wf(u_exact, u_exact, mesh,time_step, theta,theta_DG, all, DG,  false);
	CustomWeakForm wf_rhs(u_exact, u_exact, mesh,time_step, theta,theta_DG, false, DG,  true);


	CSCMatrix<double>* matrix = new CSCMatrix<double> ; 
	SimpleVector<double> * rhs = new SimpleVector<double> (ref_ndof); 
	DiscreteProblem<double> * dp = new DiscreteProblem<double> (&wf,space);	
	dp->set_linear(true,false);
	dp->assemble(matrix);
	DiscreteProblem<double> * dp_rhs = new DiscreteProblem<double> (&wf_rhs,space);	
	dp_rhs->set_linear(true,false);



  do
  {
	Hermes::Mixins::Loggable::Static::info("time=%f, ndof = %i ", current_time,ref_ndof); 
	wf_rhs.set_current_time(current_time);
	dynamic_cast<CustomInitialCondition*>(u_exact.get())->set_time(current_time);
	dp_rhs->assemble(rhs);
	UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(matrix,rhs);    
	try
	{
		solver->solve();
	}
	catch(Hermes::Exceptions::Exception e)
	{
		e.print_msg();
	}	 
     double*	vec_new = solver->get_sln_vector();
     Solution<double>::vector_to_solution(vec_new, space, u_new);

	//sview.show(u_new);
	//lview.show(u_exact);
	if(ts==1) wf_rhs.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(u_exact, u_new) );
	current_time += time_step;
	ts++;
	delete solver;
//View::wait(HERMES_WAIT_KEYPRESS);
	  }
  while (current_time<T_FINAL);


dynamic_cast<CustomInitialCondition*>(u_exact.get())->set_time(current_time);

calc_error_total(u_new, u_exact,space);


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

