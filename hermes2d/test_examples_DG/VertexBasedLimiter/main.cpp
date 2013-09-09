#include "definitions.h"
#include "../euler_util.h"
#include "algorithms.h"

const int polynomialDegree = 1;
const int initialRefinementsCount = 5;
const Algorithm algorithm = pMultigridBessiRebay;
const SolvedExample solvedExample = MovingPeak;
const EulerLimiterType limiter_type = VertexBased;

bool HermesView = true;
bool VTKView = false;

double time_step_length;
double time_interval_length;
Hermes::Mixins::Loggable logger(true);

double diffusivity = 1e-3;

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
    time_step_length = 0.01;
    time_interval_length = 1.;
    break;
  case SolidBodyRotation:
    time_step_length = 0.01;
    time_interval_length = 2 * M_PI;
    break;
  case CircularConvection:
    time_step_length = 1;
    time_interval_length = 1e4;
    break;
  case MovingPeak:
    time_step_length = 1e-3;
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
    break;
  }
  
  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition = NULL;
  ExactSolutionScalar<double>* initial_condition_der = NULL;
  ExactSolutionScalar<double>* exact_sln = NULL;
  switch(solvedExample)
  {
  case SolidBodyRotation:
    initial_condition = new InitialConditionSolidBodyRotation(mesh);
    initial_condition_der = new InitialConditionSolidBodyRotation(mesh);
    previous_initial_condition = new InitialConditionSolidBodyRotation(mesh);
    break;
  case AdvectedCube:
    initial_condition = new InitialConditionAdvectedCube(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
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
    exact_sln = new ExactSolutionMovingPeak(mesh, diffusivity, (1./2.) * M_PI);
    break;
  case Benchmark:
    initial_condition = new ZeroSolution<double>(mesh);
    initial_condition_der = new ZeroSolution<double>(mesh);
    previous_initial_condition = new ZeroSolution<double>(mesh);
    exact_sln = new ExactSolutionBenchmark(mesh, diffusivity);
    break;
  }
  
  // Solutions.
  MeshFunctionSharedPtr<double> solution(new Solution<double>);
  MeshFunctionSharedPtr<double> previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double> previous_mean_values(initial_condition);
  MeshFunctionSharedPtr<double> previous_derivatives(initial_condition_der);
  MeshFunctionSharedPtr<double> exact_solution(exact_sln);

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(520, 10, 500, 500));
  ScalarView exact_view("Exact solution", new WinGeom(520, 510, 500, 500));

  if(algorithm == MultiscaleDecomposition)
  {
    multiscale_decomposition(mesh, solvedExample, polynomialDegree, previous_mean_values, previous_derivatives, diffusivity, time_step_length,
    time_interval_length, solution, exact_solution, solution_view, exact_view);
  }
  if(algorithm == pMultigridBessiRebay)
  {
    p_multigrid(mesh, solvedExample, polynomialDegree, previous_solution, diffusivity, time_step_length, time_interval_length, 
      solution, exact_solution, solution_view, exact_view);
  }

  View::wait();
  return 0;
}
