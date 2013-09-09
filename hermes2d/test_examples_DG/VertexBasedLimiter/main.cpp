#include "definitions.h"
#include "../euler_util.h"

const int polynomialDegree = 1;
const int initialRefinementsCount = 5;
const TimeSteppingType timeSteppingType = ImplicitEuler;
const SolvedExample solvedExample = MovingPeak;
const EulerLimiterType limiter_type = VertexBased;

bool HermesView = true;
bool VTKView = false;

double time_step_length;
double time_interval_length;
Hermes::Mixins::Loggable logger(true);

double diffusivity = 1e-3;

class MyErrorCalculator : public ErrorCalculator<double>
{
public:
  MyErrorCalculator(CalculatedErrorType errorType, int component_count) : ErrorCalculator<double>(errorType)
  {
    DefaultNormFormVol<double>* form = new DefaultNormFormVol<double>(0, 0, HERMES_L2_NORM);
    this->add_error_form(form);
  }
  virtual ~MyErrorCalculator() {}
};

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
    for(int i = 1; i < al_target.cnt; i++)
      target[al_target.dof[i]] = solution_vector_fine[al_fine.dof[i-1]];
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
    time_step_length = 5e-3;
    time_interval_length = (2. * M_PI) + (time_step_length / 10.);
    break;
  case Benchmark:
    time_step_length = 1.;
    time_interval_length = 10.5 + time_step_length / 10.;
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
  case Benchmark:
    Hermes::vector<MeshSharedPtr> meshes;
    meshes.push_back(mesh);
    mloader.load("domain.msh", meshes);
    mesh->refine_all_elements();
//    mesh->refine_all_elements(2);
  //  mesh->refine_all_elements(2);
   // for(int i = 0; i < initialRefinementsCount; i++)
    //  mesh->refine_all_elements();
    break;
  }

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition;
  ExactSolutionScalar<double>* initial_condition;
  ExactSolutionScalar<double>* initial_condition_der;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
    initial_condition_der = new InitialConditionSolidBodyRotation(mesh);
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
    break;
  case AdvectedCube:
    initial_condition = new InitialConditionAdvectedCube(mesh);
    initial_condition_der = new InitialConditionAdvectedCube(mesh);
    previous_initial_condition = new InitialConditionAdvectedCube(mesh);
    break;
  case CircularConvection:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    break;
  case MovingPeak:
    initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    initial_condition_der = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    previous_initial_condition = new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.);
    break;
  case Benchmark:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    break;
  }

  MeshFunctionSharedPtr<double>previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double>previous_mean_values(initial_condition);
  MeshFunctionSharedPtr<double>previous_derivatives(initial_condition_der);

  
  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(520, 10, 500, 500));
  MeshFunctionSharedPtr<double>initial_solution(new InitialConditionBenchmark(mesh, diffusivity));
  MeshFunctionSharedPtr<double>exact_solution(new ExactSolutionMovingPeak(mesh, diffusivity, (5./2.) * M_PI));
 
  ScalarView exact_view("Exact solution", new WinGeom(520, 510, 500, 500));


#pragma region ImplicitEuler
  ImplicitWeakForm weakform_implicit(solvedExample, false, "0", "2", diffusivity);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, initial_solution));
  ExplicitWeakForm weakform_explicit(solvedExample, ExplicitEuler, 1, false, "0", "2", diffusivity);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives, initial_solution));
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
      Solution<double>::vector_to_solution(solver_explicit.get_sln_vector(), space, previous_derivatives);
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
      delete [] merged_sln;
    }
    else
    {
      solution->copy(previous_mean_values);
    }
    if(HermesView)
    {
      solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
      ((ExactSolutionMovingPeak*)(exact_solution.get()))->set_current_time(current_time + (M_PI / 2.));
      solution_view.show(solution);
      exact_view.show(exact_solution);
      MyErrorCalculator errorCalculator(RelativeErrorToGlobalNorm, 1);
      errorCalculator.calculate_errors(solution, exact_solution, true);
      std::cout << "Error: " << errorCalculator.get_error_squared(0) << std::endl;
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
