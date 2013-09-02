#include "definitions.h"
#include "../euler_util.h"

const int polynomialDegree = 1;
const int initialRefinementsCount = 7;
const TimeSteppingType timeSteppingType = ImplicitEuler;
const SolvedExample solvedExample = CircularConvection;
const EulerLimiterType limiter_type = VertexBased;

bool HermesView = true;
bool VTKView = false;

double time_step_length;
double time_interval_length;
Hermes::Mixins::Loggable logger(true);

//int main(int argc, char* argv[])
//{
//  return test();
//}

int main(int argc, char* argv[])
{
  // test();
  Hermes::Mixins::Loggable::set_static_logFile_name("logfile.h2d");

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
    time_step_length = 1e-1;
    time_interval_length = 1e4;
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
  }

  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  // Previous time level solution (initialized by the initial condition).
  ExactSolutionScalar<double>* previous_initial_condition;
  ExactSolutionScalar<double>* initial_condition;
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
    previous_initial_condition = new ZeroSolution<double>(mesh);
    break;
  }

  MeshFunctionSharedPtr<double>previous_solution(previous_initial_condition);
  MeshFunctionSharedPtr<double>previous_solution_time_step(initial_condition);
  MeshFunctionSharedPtr<double>updated_previous_solution_time_step(initial_condition);
  MeshFunctionSharedPtr<double>exact_solution_circular(new InitialConditionCircularConvection(mesh));

  // Visualization.
  ScalarView solution_view("Solution", new WinGeom(520, 10, 500, 500));

#pragma region ExplicitRK
  // Weak form.
  ExplicitWeakForm weakform_1(solvedExample, timeSteppingType, 1, solvedExample == CircularConvection);
  ExplicitWeakForm weakform_2(solvedExample, timeSteppingType, 2, solvedExample == CircularConvection);
  ExplicitWeakForm weakform_3(solvedExample, timeSteppingType, 3, solvedExample == CircularConvection);
  weakform_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution_time_step, exact_solution_circular));
  weakform_1.set_current_time_step(time_step_length);
  weakform_2.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution, previous_solution_time_step, exact_solution_circular));
  weakform_2.set_current_time_step(time_step_length);
  weakform_3.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution, previous_solution_time_step, exact_solution_circular));
  weakform_3.set_current_time_step(time_step_length);

  // Solver.
  LinearSolver<double> solver_1(&weakform_1, space);
  LinearSolver<double> solver_2(&weakform_2, space);
  LinearSolver<double> solver_3(&weakform_3, space);
  solver_1.set_jacobian_constant();
  solver_2.set_jacobian_constant();
  solver_3.set_jacobian_constant();
#pragma endregion

#pragma region ImplicitEuler
  ImplicitWeakForm weakform_implicit(solvedExample, true, "Inlet", "Bdy");
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_solution_time_step, exact_solution_circular));
  ExplicitWeakForm weakform_explicit(solvedExample, ExplicitEuler, 1, solvedExample == CircularConvection, "Inlet", "Bdy");
  weakform_explicit.set_current_time_step(time_step_length);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(updated_previous_solution_time_step, exact_solution_circular));
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
    if(timeSteppingType == ExplicitRK)
    {
      // 1st step.
      solver_1.solve();
      PostProcessing::Limiter<double>* limiter_1 = create_limiter(limiter_type, space, solver_1.get_sln_vector(), polynomialDegree);
      solution->copy(limiter_1->get_solution());

      // 2nd step.
      solver_2.solve();
      PostProcessing::Limiter<double>* limiter_2 = create_limiter(limiter_type, space, solver_2.get_sln_vector(), polynomialDegree);
      previous_solution->copy(limiter_2->get_solution());

      // 3rd step.
      solver_3.solve();
      PostProcessing::Limiter<double>* limiter_3 = create_limiter(limiter_type, space, solver_3.get_sln_vector(), polynomialDegree);
      solution->copy(limiter_3->get_solution());
    }
    else if(timeSteppingType == ImplicitEuler)
    {
      if(time_step)
      {
        // 0th - step
        double* previous_sln_vector = new double[ndofs];
        double* sln_vector;
        double* mean_values;
        OGProjection<double>::project_global(space, previous_solution_time_step, previous_sln_vector);
        // 1 - solve implicit.
        solver_implicit.solve();
        mean_values = new double[ndofs];
        Solution<double>::vector_to_solution(solver_implicit.get_sln_vector(), const_space, updated_previous_solution_time_step);
        OGProjection<double>::project_global(space, updated_previous_solution_time_step, mean_values);

        // 2 - Update the mean values.
        Element* e;
        for_all_active_elements(e, mesh)
        {
          AsmList<double> al_fine;
          space->get_element_assembly_list(e, &al_fine);
          for(unsigned int shape_i = 0; shape_i < al_fine.cnt; shape_i++)
          {
            int order = space->get_shapeset()->get_order(al_fine.idx[shape_i], e->get_mode());
            if(order == 0)
            {
              int dof = al_fine.dof[shape_i];
              previous_sln_vector[dof] = mean_values[dof];
            }
          }
        }

        // Clean up.
        delete [] mean_values;

        // 3 - Solve explicit.
        Solution<double>::vector_to_solution(previous_sln_vector, space, updated_previous_solution_time_step);
        delete [] previous_sln_vector;
      }
      solver_explicit.solve();
      if(polynomialDegree)
      {
        PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, space, solver_explicit.get_sln_vector(), polynomialDegree);
        solution->copy(limiter->get_solution());
        delete limiter;
      }
    }

    if(HermesView)
    {
      solution_view.set_title("Solution - time step: %i, time: %f.", time_step, current_time);
      solution_view.show(solution);
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

    previous_solution_time_step->copy(solution);
    current_time += time_step_length;
  }

  View::wait();
  return 0;
}