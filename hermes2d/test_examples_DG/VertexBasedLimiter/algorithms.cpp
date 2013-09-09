#include "algorithms.h"

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_mean_values, 
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double time_step_length, 
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                              ScalarView solution_view, ScalarView exact_view)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));

  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs = space->get_num_dofs();

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  ImplicitWeakForm weakform_implicit(solvedExample, false, "0", "2", diffusivity);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_implicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives));
  ExplicitWeakForm weakform_explicit(solvedExample, false, "0", "2", diffusivity);
  weakform_explicit.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_mean_values, previous_derivatives));
  weakform_explicit.set_current_time_step(time_step_length);
  LinearSolver<double> solver_implicit(&weakform_implicit, const_space);
  LinearSolver<double> solver_explicit(&weakform_explicit, space);

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

    solution_view.show(solution);
    
    /*
    ((ExactSolutionMovingPeak*)(exact_solution.get()))->set_current_time(current_time + (M_PI / 2.));
    exact_view.show(exact_solution);
    MyErrorCalculator errorCalculator(RelativeErrorToGlobalNorm, 1);
    errorCalculator.calculate_errors(solution, exact_solution, true);
    std::cout << "Error: " << errorCalculator.get_error_squared(0) << std::endl;
    */

    /*
    {
    Linearizer lin;
    char* filename = new char[100];
    sprintf(filename, "Sln-%i.vtk", time_step);
    lin.save_solution_vtk(solution, filename, "sln");
    delete [] filename;
    }
    */
    current_time += time_step_length;
  }
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                 double diffusivity, double time_step_length, 
                 double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution, 
                 ScalarView solution_view, ScalarView exact_view)
{
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);

  // Spaces
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Previous iteration solution
  MeshFunctionSharedPtr<double> prev_iter_solution(new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.));

  // 1 - solver
  SmoothingWeakForm weakform_1(solvedExample);
  weakform_1.set_current_time_step(time_step_length);
  weakform_1.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_iter_solution, previous_sln));
  LinearSolver<double> solver_1(&weakform_1, space_1);

  // 1 - Residual measurement.
  MeshFunctionSharedPtr<double> aux_solution(new ExactSolutionMovingPeak(mesh, diffusivity, M_PI / 2.));
  SmoothingWeakFormResidual weakform_residual(solvedExample);
  weakform_residual.set_current_time_step(time_step_length);
  weakform_residual.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_iter_solution, aux_solution));
  DiscreteProblem<double> dp(&weakform_residual, space_1);
  Algebra::UMFPackVector<double> vec;

  // 0 - solver
  CoarseWeakForm weakform_0(solvedExample);
  weakform_0.set_current_time_step(time_step_length);
  weakform_0.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_iter_solution, previous_sln));
  LinearSolver<double> solver_0(&weakform_0, space_0);

  double current_time = 0.;
  int number_of_steps = (time_interval_length - current_time) / time_step_length;
  for(int time_step = 0; time_step <= number_of_steps; time_step++)
  { 
    std::cout << "Time step " << time_step << std::endl;

    // 1st part - smoothing on the 1st level.
    dp.set_space(space_1);
    for(int iteration_1 = 1; iteration_1 < 15; iteration_1++)
    {
      solver_1.solve();
      Solution<double>::vector_to_solution(solver_1.get_sln_vector(), space_1, prev_iter_solution);
      solution_view.show(prev_iter_solution);
      
      dp.assemble(&vec);
      aux_solution->copy(prev_iter_solution);
      double residual_norm = Hermes2D::get_l2_norm(&vec);
      std::cout << "\tIteration - 1: " << iteration_1 << ", residual norm: " << residual_norm << std::endl;
      if(residual_norm < 1e-4)
        break;
    }

    // Take only the 0 - part of the previous solution
    aux_solution->copy(prev_iter_solution);
    dp.set_space(space_0);
    for(int iteration_0 = 1; iteration_0 < 15; iteration_0++)
    {
      solver_0.solve();
      Solution<double>::vector_to_solution(solver_0.get_sln_vector(), space_0, prev_iter_solution);
      solution_view.show(prev_iter_solution);
      
      dp.assemble(&vec);
      aux_solution->copy(prev_iter_solution);
      double residual_norm = Hermes2D::get_l2_norm(&vec);
      std::cout << "\tIteration - 2: " << iteration_0 << ", residual norm: " << residual_norm << std::endl;
      if(residual_norm < 1e-4)
        break;
    }

    double* solution_vector = merge_slns(solver_0.get_sln_vector(), space_0, solver_1.get_sln_vector(), space_1, space_1);
    Solution<double>::vector_to_solution(solution_vector, space_1, previous_sln);
    
    solution_view.show(previous_sln);

    ((ExactSolutionMovingPeak*)(exact_solution.get()))->set_current_time(current_time + (M_PI / 2.));
    exact_view.show(exact_solution);
    MyErrorCalculator errorCalculator(RelativeErrorToGlobalNorm, 1);
    errorCalculator.calculate_errors(previous_sln, exact_solution, true);
    std::cout << "Error: " << errorCalculator.get_error_squared(0) << std::endl;

    /*
    {
    Linearizer lin;
    char* filename = new char[100];
    sprintf(filename, "Sln-%i.vtk", time_step);
    lin.save_solution_vtk(solution, filename, "sln");
    delete [] filename;
    }
    */
    current_time += time_step_length;
  }
}