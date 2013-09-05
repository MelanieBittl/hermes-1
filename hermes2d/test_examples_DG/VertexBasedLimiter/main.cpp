#include "definitions.h"
#include "../euler_util.h"

const int polynomialDegree = 1;
const int initialRefinementsCount = 7;
const TimeSteppingType timeSteppingType = ImplicitEuler;
const SolvedExample solvedExample = MovingPeak;
const EulerLimiterType limiter_type = VertexBased;

bool HermesView = true;
bool VTKView = true;

double time_step_length;
double time_interval_length;
Hermes::Mixins::Loggable logger(true);

double diffusivity = 1e-4;

//int main(int argc, char* argv[])
//{
//  return test();
//}

double* merge_slns(double* solution_vector_coarse, SpaceSharedPtr<double> space_coarse, double* solution_vector_fine, SpaceSharedPtr<double> space_fine, SpaceSharedPtr<double> space_full)
{
  double* target = new double[space_full->get_num_dofs()];
  Element *e;
  for_all_active_elements(e, space_full->get_mesh())
  {
    AsmList<double> al_coarse, al_fine, al_target;
    space_coarse->get_element_assembly_list(e, &al_coarse);
    space_fine->get_element_assembly_list(e, &al_fine);
    space_full->get_element_assembly_list(e, &al_target);

    target[al_target.dof[0]] = solution_vector_coarse[al_coarse.dof[0]];
    target[al_target.dof[1]] = solution_vector_fine[al_fine.dof[0]];
    target[al_target.dof[2]] = solution_vector_fine[al_fine.dof[1]];
  }
  return target;
}

int main(int argc, char* argv[])
{
  if(argc > 1)
    diffusivity = atof(argv[1]);
  // test();
  Hermes::Mixins::Loggable::set_static_logFile_name("logfile.h2d");
  HermesCommonApi.set_integral_param_value(numThreads, 16);

  switch(solvedExample)
  {
  case AdvectedCube:
    time_step_length = (timeSteppingType == ExplicitRK || timeSteppingType == ExplicitEuler) ? 0.001 : 0.01;
    time_interval_length = 1.;
    break;
  case SolidBodyRotation:
    time_step_length = (timeSteppingType == ExplicitRK || timeSteppingType == ExplicitEuler) ? 0.001 : 0.01;
    time_interval_length = 2 * M_PI;
    break;
  case CircularConvection:
    time_step_length = 1;
    time_interval_length = 1e4;
    break;
  case MovingPeak:
    time_step_length = 1e-1;
    time_interval_length = 5. * M_PI / 2.;
    break;
  }

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    mloader.load("domain_rotation.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case AdvectedCube:
    mloader.load("larger_domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case CircularConvection:
    mloader.load("domain_circular_convection.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  case MovingPeak:
    mloader.load("domain.xml", mesh);
    for(int i = 0; i < initialRefinementsCount; i++)
      mesh->refine_all_elements();
    break;
  }

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition;
  ExactSolutionScalar<double>* updated_previous_initial_condition;
  ExactSolutionScalar<double>* initial_condition;
  ExactSolutionScalar<double>* initial_condition_der;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
    break;
  case AdvectedCube:
    initial_condition = new InitialConditionAdvectedCube(mesh);
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);
    break;
  case CircularConvection:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    updated_previous_initial_condition = new ZeroSolution<double>(mesh);
    break;
    case MovingPeak:
    initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    initial_condition_der = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    previous_initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    updated_previous_initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    break;
  }

  MeshFunctionSharedPtr<double>previous_solution(previous_initial_condition);

  MeshFunctionSharedPtr<double>previous_mean_values(initial_condition);

  MeshFunctionSharedPtr<double>previous_derivatives(initial_condition_der);
  MeshFunctionSharedPtr<double>updated_previous_mean_values(updated_previous_initial_condition);
  MeshFunctionSharedPtr<double>exact_solution_circular(new InitialConditionCircularConvection(mesh));

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(520, 10, 500, 500));

#pragma region ImplicitEuler
  ImplicitWeakForm weakform_implicit(solvedExample, false, "Inlet", "Bdy", diffusivity);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, exact_solution_circular));
  ExplicitWeakForm weakform_explicit(solvedExample, ExplicitEuler, 1, false, "Inlet", "Bdy", diffusivity);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, exact_solution_circular));
  weakform_explicit.set_current_time_step(time_step_length);
  LinearSolver<double> solver_implicit(&weakform_implicit, const_space);
  LinearSolver<double> solver_explicit(&weakform_explicit, space);
#pragma endregion

  // Solution.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);

  ScalarView sol_view;

  double current_time = 0.;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
  for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    std::cout << "Iteration " << time_step << std::endl;
    solver_implicit.solve();
    Solution<double>::vector_to_solution(solver_implicit.get_sln_vector(), const_space, previous_mean_values);

    if(polynomialDegree)
    {
      solver_explicit.solve();
      double* merged_sln = merge_slns(solver_implicit.get_sln_vector(), const_space,solver_explicit.get_sln_vector(), space, full_space);

      PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, full_space, merged_sln, polynomialDegree);
      solution->copy(limiter->get_solution());
      for(int i = 0; i < full_space->get_num_dofs(); i++)
        if(!(i%3))
          limiter->get_solution_vector()[i] = 0;
      Solution<double>::vector_to_solution(limiter->get_solution_vector(), full_space, previous_derivatives);
      delete limiter;
    
      delete [] merged_sln;
    }
    else
    {
      solution->copy(previous_mean_values);
    }
    if(HermesView)
    {
      solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
      solution_view.show(solution);
      DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
      errorCalculator.calculate_errors(solution, exact_solution_circular, true);
      std::cout << "Error: " << errorCalculator.get_error_squared(0) << std::endl;

      //View::wait_for_keypress();
    }
    if(VTKView)
    {
      Linearizer lin;
      char* filename = new char[100];
      sprintf(filename, "Sln-%i.vtk", time_step);
      lin.save_solution_vtk(solution, filename, "sln");
      delete [] filename;
    }

    current_time += time_step_length;
  }

  View::wait();
  return 0;
}
