#define HERMES_REPORT_ALL
#include "definitions.h"
#include "hermes2d.h"
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.



MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";


int main(int argc, char* argv[])
{ 
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements(); 
Hermes2DApi.set_integral_param_value(numThreads, 1); 
  //L2Space
//L2_SEMI_CG_Space
  // Create an L2 	space with default shapeset.
  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh,P_INIT));	
/*BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
bview.show(space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));
 
  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView s2view("solution_2", new WinGeom(0, 500, 500, 400));
s2view.set_min_max_range(0, 0.5);
	ScalarView fview("filter", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	lview.show(u_prev_time);

	
	int ndof = space->get_num_dofs();


CustomWeakForm wf(u_prev_time,true, true);
 // Initialize linear solver.
 Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, u_new);
		s2view.show(u_new);
  }catch(std::exception& e)
  {
    std::cout << e.what();
  }
 






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

