#include "definitions.h"

const int polynomialDegree = 2;
const int initialRefinementsCount = 5;
const double time_step_length = 0.015;
const double time_interval_length = 1.;
const double logPercentTimeSteps = 1.;
int logPeriod = (int)std::max<double>(1., ((logPercentTimeSteps / 100.) * (time_interval_length / time_step_length)));

Hermes::Mixins::Loggable logger(true);

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("larger_domain.xml", mesh);

  for(int i = 0; i < initialRefinementsCount; i++)
    mesh->refine_all_elements();

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  // Visualization.
  // BaseView<double> space_view("Space browser", new WinGeom(10, 10, 500, 500));
  // space_view.show(space, HERMES_EPS_HIGH);
  // View::wait_for_keypress();

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>initial_condition(new CustomInitialCondition(mesh));
  // Visualization.
  ScalarView solution_view("Initial condition", new WinGeom(520, 10, 500, 500));
  solution_view.show(initial_condition);

  // Weak form.
  CustomWeakForm weakform(true);
  weakform.set_ext(initial_condition);
  weakform.set_current_time_step(time_step_length);

  // Solver.
  LinearSolver<double> solver(&weakform, space);

  // Solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);

  double current_time = 0.;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
  for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    if((!(time_step % logPeriod)) || (time_step == number_of_steps))
    {
      logger.info("Time step: %i, time: %f.", time_step, current_time);
      solver.set_verbose_output(true);
    }

    solver.solve();
    //Solution<double>::vector_to_solution(solver.get_sln_vector(), space, solution);
    PostProcessing::VertexBasedLimiter limiter(space, solver.get_sln_vector(), polynomialDegree);
    solution = limiter.get_solution();
    
    if((!(time_step % logPeriod)) || (time_step == number_of_steps))
    {
      solver.set_verbose_output(false);
      solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
      solution_view.show(solution, HERMES_EPS_NORMAL, H2D_FN_DX_0);
      //View::wait_for_keypress();
    }

    initial_condition->copy(solution);
    current_time += time_step_length;
  }

  View::wait();
  return 0;
}