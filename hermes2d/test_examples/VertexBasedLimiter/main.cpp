#include "definitions.h"

const int polynomialDegree = 2;
const int initialRefinementsCount = 6;
const double logPercentTimeSteps = 1.;
const TimeSteppingType timeSteppingType = Explicit;
const SolvedExample solvedExample = SolidBodyRotation;

bool HermesView = false;
bool VTKView = true;

const double time_step_length = timeSteppingType == Explicit ? 0.001 : 0.01;
const double time_interval_length = solvedExample == SolidBodyRotation ? 2 * M_PI : 1.;
int logPeriod = (int)std::max<double>(1., ((logPercentTimeSteps / 100.) * (time_interval_length / time_step_length)));
Hermes::Mixins::Loggable logger(true);

//int main(int argc, char* argv[])
//{
//  return test();
//}

int main(int argc, char* argv[])
{
  // test();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load((solvedExample == SolidBodyRotation ? "domain_rotation.xml" : "larger_domain.xml"), mesh);

  for(int i = 0; i < initialRefinementsCount; i++)
    mesh->refine_all_elements();

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* initial_condition;
  if(solvedExample == SolidBodyRotation)
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
  else if (solvedExample == AdvectedCube)
    initial_condition = new InitialConditionAdvectedCube(mesh);

  ExactSolutionScalar<double>* previous_initial_condition;
  if(solvedExample == SolidBodyRotation)
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
  else if (solvedExample == AdvectedCube)
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);

  MeshFunctionSharedPtr<double>previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double>previous_solution_time_step(initial_condition);
  // Visualization.
  ScalarView solution_view("Initial condition", new WinGeom(520, 10, 500, 500));
  Linearizer lin;
  if(HermesView)
    solution_view.show(previous_solution);

  // Weak form.
  CustomWeakForm weakform_1(solvedExample, timeSteppingType, 1);
  CustomWeakForm weakform_2(solvedExample, timeSteppingType, 2);
  CustomWeakForm weakform_3(solvedExample, timeSteppingType, 3);
  weakform_1.set_ext(previous_solution_time_step);
  weakform_1.set_current_time_step(time_step_length);
  
  weakform_2.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution, previous_solution_time_step));
  weakform_2.set_current_time_step(time_step_length);
  
  weakform_3.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution, previous_solution_time_step));
  weakform_3.set_current_time_step(time_step_length);

  // Solver.
  LinearSolver<double> solver_1(&weakform_1, space);
  LinearSolver<double> solver_2(&weakform_2, space);
  LinearSolver<double> solver_3(&weakform_3, space);
  if(timeSteppingType == Explicit)
  {
    solver_1.set_jacobian_constant();
    solver_2.set_jacobian_constant();
    solver_3.set_jacobian_constant();
  }

  Hermes::Mixins::Loggable::set_static_logFile_name("logfile.h2d");

  // Solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);

  double current_time = 0.;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
	for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    if((!(time_step % logPeriod)) || (time_step == number_of_steps))
    {
      logger.info("Time step: %i, time: %f.", time_step, current_time);
      solver_1.set_verbose_output(true);
      solver_2.set_verbose_output(true);
      solver_3.set_verbose_output(true);
    }

    // 1st step.
    solver_1.solve();
    PostProcessing::VertexBasedLimiter limiter_1(space, solver_1.get_sln_vector(), polynomialDegree);
    previous_solution->copy(limiter_1.get_solution());

    // 2nd step.
    solver_2.solve();
    PostProcessing::VertexBasedLimiter limiter_2(space, solver_2.get_sln_vector(), polynomialDegree);
    previous_solution->copy(limiter_2.get_solution());

    // 3rd step.
    solver_3.solve();
    PostProcessing::VertexBasedLimiter limiter_3(space, solver_3.get_sln_vector(), polynomialDegree);
    solution->copy(limiter_3.get_solution());

    if((!(time_step % logPeriod)) || (time_step == number_of_steps))
    {
      solver_1.set_verbose_output(false);
      solver_2.set_verbose_output(false);
      solver_3.set_verbose_output(false);
      if(HermesView)
      {
        solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
        solution_view.show(solution);
        //View::wait_for_keypress();
      }
      if(VTKView)
      {
        char* filename = new char[100];
        sprintf(filename, "Sln-%i.vtk", time_step);
        lin.save_solution_vtk(solution, filename, "sln");
        delete [] filename;
      }
    }

    previous_solution_time_step->copy(solution);
    current_time += time_step_length;
  }

  View::wait();
  return 0;
}