#include "algorithms.h"

static double exact_solver_error;
static const double tolerance = 1e-4;
double initial_error = -1;
MeshFunctionSharedPtr<double> es(new Solution<double>());

double calc_l2_error(SolvedExample solvedExample, MeshSharedPtr mesh, MeshFunctionSharedPtr<double> fn_1, MeshFunctionSharedPtr<double> fn_2, Hermes::Mixins::Loggable& logger)
{
  ErrorWeakForm wf(solvedExample);
  SpaceSharedPtr<double> mspace(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  wf.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(fn_1, fn_2));
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, mspace);
  SimpleVector<double> vector;
  dp->assemble(&vector);
  double result = 0.;
  for(int i = 0; i < vector.get_size(); i++)
    result += vector.get(i);
  result = std::sqrt(result);
  logger.info("L2 Error: %g.", result);
  return result;
}

bool error_condition(double error)
{
  return std::abs(error / initial_error) < tolerance;
}

void solve_exact(SolvedExample solvedExample, SpaceSharedPtr<double> space, double diffusivity, double s, double sigma, MeshFunctionSharedPtr<double> exact_solution, MeshFunctionSharedPtr<double> initial_sln, double time_step, Hermes::Mixins::Loggable& logger)
{
  MeshFunctionSharedPtr<double> exact_solver_sln(new Solution<double>());
  ScalarView* exact_solver_view = new ScalarView("Exact solver solution", new WinGeom(0, 360, 600, 350));

  // Exact solver
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", diffusivity, s, sigma, initial_sln);
  weakform_exact.set_current_time_step(time_step);
  LinearSolver<double> solver_exact(&weakform_exact, space);
  solver_exact.solve();
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, exact_solver_sln);
  exact_solver_error = calc_l2_error(solvedExample, space->get_mesh(), exact_solver_sln, exact_solution, logger);
  exact_solver_view->show(exact_solver_sln);
  initial_error = get_l2_norm(solver_exact.get_sln_vector(), space->get_num_dofs());
  Solution<double>::vector_to_solution(solver_exact.get_sln_vector(), space, es);
}

void multiscale_decomposition(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_mean_values,
                              MeshFunctionSharedPtr<double> previous_derivatives, double diffusivity, double s, double sigma, double time_step_length,
                              double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
                              ScalarView* solution_view, ScalarView* exact_view, Hermes::Mixins::Loggable& logger)
{
  // Standard L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor(false)));
  SpaceSharedPtr<double> full_space(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));

  int ndofs = space->get_num_dofs();
  int const_ndofs = const_space->get_num_dofs();
  int full_ndofs = full_space->get_num_dofs();

  OGProjection<double>::project_global(const_space, previous_mean_values, previous_mean_values);
  OGProjection<double>::project_global(space, previous_derivatives, previous_derivatives);

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", diffusivity, s, sigma, exact_solution);
  MultiscaleWeakForm weakform_implicit(solvedExample, true, "Inlet", diffusivity, s, sigma, exact_solution, false);
  MultiscaleWeakForm weakform_explicit(solvedExample, true, "Inlet", diffusivity, s, sigma, exact_solution, true);
  ExplicitWeakForm weakform_explicit_offdiag(solvedExample, true, "Inlet", diffusivity, s, sigma);
  MassWeakForm weakform_mass;
  weakform_exact.set_current_time_step(time_step_length);
  weakform_implicit.set_current_time_step(time_step_length);
  weakform_explicit.set_current_time_step(time_step_length);
  weakform_explicit_offdiag.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_full;
  CSCMatrix<double> matrix_A_der;
  SimpleVector<double> vector_b_der;
  CSCMatrix<double> matrix_M_der;
  CSCMatrix<double> matrix_A_offdiag;

  CSCMatrix<double> matrix_A_means;
  SimpleVector<double> vector_b_means;
  CSCMatrix<double> matrix_M_means;

  // Assembler.
  DiscreteProblem<double> dp;
  // Level 2.
  dp.set_space(full_space);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_full);

  // Level 1.
  dp.set_space(space);
  dp.set_weak_formulation(&weakform_explicit);
  dp.assemble(&matrix_A_der, &vector_b_der);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_der);
  dp.set_weak_formulation(&weakform_explicit_offdiag);
  dp.assemble(&matrix_A_offdiag);

  // Level 0.
  dp.set_space(const_space);
  dp.set_weak_formulation(&weakform_implicit);
  dp.assemble(&matrix_A_means, &vector_b_means);
  dp.set_weak_formulation(&weakform_mass);
  dp.assemble(&matrix_M_means);

  SimpleVector<double> vector_A_der(ndofs);
  SimpleVector<double> vector_A_means(const_ndofs);

  UMFPackLinearMatrixSolver<double> solver_means(&matrix_A_means, &vector_A_means);
  solver_means.setup_factorization();
  solver_means.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_der(&matrix_A_der, &vector_A_der);
  solver_der.setup_factorization();
  solver_der.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  SimpleVector<double> sln_means(const_ndofs);
  SimpleVector<double> sln_means_long(full_ndofs);
  SimpleVector<double> sln_means_long_temp(full_ndofs);
  SimpleVector<double> sln_der(ndofs);
  SimpleVector<double> sln_der_long(full_ndofs);
  SimpleVector<double> sln_der_long_temp(full_ndofs);
  SimpleVector<double> sln_der_offdiag(ndofs);

  OGProjection<double>::project_global(const_space, previous_mean_values, sln_means.v);
  OGProjection<double>::project_global(space, previous_derivatives, sln_der.v);

  // Reporting.
  int num_coarse = 0;
  int num_fine = 0;
  int iterations = 0;

  for(int iteration = 1;; iteration++)
  {
    iterations++;
    logger.info("Iteration %i.", iteration);
    num_coarse++;

    matrix_M_means.multiply_with_vector(sln_means.v, vector_A_means.v, true);
    add_means(&sln_der, &sln_der_long, space, full_space);
    matrix_A_full.multiply_with_vector(sln_der_long.v, sln_der_long_temp.v, true);
    vector_A_means.add_vector(cut_off_ders(sln_der_long_temp.v, const_space, full_space)->change_sign());
    vector_A_means.add_vector(&vector_b_means);
    solver_means.solve();
    sln_means.set_vector(solver_means.get_sln_vector());
    //Solution<double>::vector_to_solution(solver_means.get_sln_vector(), const_space, previous_mean_values);
    //solution_view->show(previous_mean_values);
    //solution_view->wait_for_keypress();

    if(polynomialDegree)
    {
      num_fine++;

      matrix_M_der.multiply_with_vector(sln_der.v, vector_A_der.v, true);

      add_ders(&sln_means, &sln_means_long, const_space, full_space);
      matrix_A_full.multiply_with_vector(sln_means_long.v, sln_means_long_temp.v, true);
      vector_A_der.add_vector(cut_off_means(sln_means_long_temp.v, space, full_space)->change_sign());

      matrix_A_offdiag.multiply_with_vector(sln_der.v, sln_der_offdiag.v, true);
      vector_A_der.add_vector(sln_der_offdiag.change_sign());

      vector_A_der.add_vector(&vector_b_der);
      solver_der.solve();
      sln_der.set_vector(solver_der.get_sln_vector());

      double* merged_sln = merge_slns(sln_means.v, const_space, sln_der.v, space, full_space);
      //Solution<double>::vector_to_solution(sln_der.v, space, previous_derivatives);
      //solution_view->show(previous_derivatives);
      //solution_view->wait_for_keypress();
      Solution<double>::vector_to_solution(merged_sln, full_space, solution);
      delete [] merged_sln;
    }
    else
    {
      solution->copy(previous_mean_values);
    }

    solution_view->show(solution);
    //solution_view->wait_for_keypress();

    if(error_condition(calc_l2_error(solvedExample, mesh, solution, es, logger)))
      break;
  }

  logger.info("Iterations: %i", iterations);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("Fine systems solved: %i", num_fine);
}

bool show_intermediate = false;

static double residual_drop(int level, bool pre)
{
  return 1e-1;
}

bool residual_condition(CSCMatrix<double>* mat, SimpleVector<double>* vec, double* sln_vector, double* residual, Hermes::Mixins::Loggable& logger, int iteration, bool output)
{
  mat->multiply_with_vector(sln_vector, residual, true);
  for(int i = 0; i < mat->get_size(); i++)
    residual[i] -= vec->get(i);

  if(output)
  {
    double residual_norm = Hermes2D::get_l2_norm(residual, mat->get_size());
    logger.info("\tIteration: %i, residual norm: %g.", iteration, residual_norm);
  }

  return false;
}

void p_multigrid(MeshSharedPtr mesh, SolvedExample solvedExample, int polynomialDegree, MeshFunctionSharedPtr<double> previous_sln,
                 double diffusivity, double time_step_length,
                 double time_interval_length, MeshFunctionSharedPtr<double> solution, MeshFunctionSharedPtr<double> exact_solution,
                 ScalarView* solution_view, ScalarView* exact_view, double s, double sigma, Hermes::Mixins::Loggable& logger, int steps)
{
  // Spaces
  SpaceSharedPtr<double> space_2(new L2Space<double>(mesh, polynomialDegree, new L2ShapesetTaylor));
  int ndofs_2 = space_2->get_num_dofs();
  SpaceSharedPtr<double> space_1(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  int ndofs_1 = space_1->get_num_dofs();
  SpaceSharedPtr<double> space_0(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  int ndofs_0 = space_0->get_num_dofs();

  // Matrices A, vectors b.
  ExactWeakForm weakform_exact(solvedExample, true, "Inlet", diffusivity, s, sigma, exact_solution);
  weakform_exact.set_current_time_step(time_step_length);
  CSCMatrix<double> matrix_A_2;
  SimpleVector<double> vector_b_2;
  CSCMatrix<double> matrix_A_1;
  SimpleVector<double> vector_b_1;
  CSCMatrix<double> matrix_A_0;
  SimpleVector<double> vector_b_0;

  // Matrices (M+A_tilde), vectors -A(u_K)
  SmoothingWeakForm weakform_smoother(solvedExample, true, 1, true, "Inlet", diffusivity, s, sigma);
  SmoothingWeakForm weakform_smoother_coarse(solvedExample, false, 1, true, "Inlet", diffusivity, s, sigma);
  weakform_smoother.set_current_time_step(time_step_length);
  weakform_smoother.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  weakform_smoother_coarse.set_current_time_step(time_step_length);
  weakform_smoother_coarse.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(previous_sln, exact_solution));
  CSCMatrix<double> matrix_MA_tilde_2;
  SimpleVector<double> vector_A_2(ndofs_2);
  CSCMatrix<double> matrix_MA_tilde_1;
  SimpleVector<double> vector_A_1(ndofs_1);
  CSCMatrix<double> matrix_MA_0;
  SimpleVector<double> vector_A_0(ndofs_0);

  // Assembler.
  DiscreteProblem<double> dp;
  // Level 2.
  dp.set_space(space_2);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_2, &vector_b_2);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_2);

  // Level 1.
  dp.set_space(space_1);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_1, &vector_b_1);
  dp.set_weak_formulation(&weakform_smoother);
  dp.assemble(&matrix_MA_tilde_1);

  // Level 0.
  dp.set_space(space_0);
  dp.set_weak_formulation(&weakform_exact);
  dp.assemble(&matrix_A_0, &vector_b_0);
  dp.set_weak_formulation(&weakform_smoother_coarse);
  dp.assemble(&matrix_MA_0);

  UMFPackLinearMatrixSolver<double> solver_2(&matrix_MA_tilde_2, &vector_A_2);
  solver_2.setup_factorization();
  solver_2.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_1(&matrix_MA_tilde_1, &vector_A_1);
  solver_1.setup_factorization();
  solver_1.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);
  UMFPackLinearMatrixSolver<double> solver_0(&matrix_MA_0, &vector_A_0);
  solver_0.setup_factorization();
  solver_0.set_reuse_scheme(HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY);

  // Utils.
  double* residual_2 = new double[ndofs_2];
  SimpleVector<double> sln_2(ndofs_2);
  double* residual_1 = new double[ndofs_1];
  SimpleVector<double> sln_1(ndofs_1);
  double* residual_0 = new double[ndofs_0];
  SimpleVector<double> sln_0(ndofs_0);

  double initial_residual_norm;
  // Reports.
  int num_coarse = 0;
  int num_2 = 0;
  int num_1 = 0;
  int v_cycles = 0;

  for(int step = 1;; step++)
  {
    logger.info("V-cycle %i.", step);
    v_cycles++;

#pragma region 0 - highest level
    // Store the previous solution.
    OGProjection<double>::project_global(space_2, previous_sln, &sln_2);
    if(polynomialDegree > 1)
    {
      for(int iteration = 1; iteration <= steps; iteration++)
      {
        // Solve for increment.
        matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
        vector_A_2.change_sign()->add_vector(&vector_b_2);
        solver_2.solve();
        sln_2.add_vector(solver_2.get_sln_vector());

        // Make solution
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->show(previous_sln);

        // Residual check.
        residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger, iteration, true);
      }
    }
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    OGProjection<double>::project_global(space_1, previous_sln, &sln_1);

    // f_P1
    SimpleVector<double> f_P1(ndofs_1);
    f_P1.zero();
    // Minus A_P1
    SimpleVector<double> R_P1(ndofs_1);
    // Minus(minus) projected_A_P1
    SimpleVector<double> projected_A_2(ndofs_2);
    matrix_A_2.multiply_with_vector(sln_2.v, projected_A_2.v, true);

    SimpleVector<double>* projected_A_P_1
      = (SimpleVector<double>*)cut_off_quadratic_part(projected_A_2.v, space_1, space_2);

    SimpleVector<double>* sln_2_projected = cut_off_quadratic_part(sln_2.v, space_1, space_2);
    matrix_A_1.multiply_with_vector(sln_2_projected->v, R_P1.v, true);

    R_P1.change_sign();
    f_P1.add_vector(&R_P1);
    f_P1.add_vector(projected_A_P_1);
    f_P1.change_sign();

    for(int iteration = 1; iteration <= steps; iteration++)
    {
      // Solve for increment.
      SimpleVector<double>* rhs;
      if(iteration == 1)
      {
        if(polynomialDegree > 1)
        {
          rhs = cut_off_quadratic_part(residual_2, space_1, space_2);
          memcpy(vector_A_1.v, rhs->v, ndofs_1 * sizeof(double));
        }
        else
        {
          matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);
          vector_A_1.change_sign()->add_vector(&vector_b_1);
        }
      }
      else
      {
        // A(u_K) - done after the first step.
        matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

        if(polynomialDegree > 1)
          vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
        else
          vector_A_1.change_sign()->add_vector(&vector_b_1);
      }
      solver_1.solve();
      sln_1.add_vector(solver_1.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      if(show_intermediate)
        solution_view->show(previous_sln);

      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, false);
    }
#pragma endregion

#pragma region 2 - Solve the problem on the coarse level exactly
    OGProjection<double>::project_global(space_0, previous_sln, &sln_0);

    // f_P0
    SimpleVector<double> f_P0(ndofs_0);
    f_P0.zero();
    // Minus A_P0
    SimpleVector<double> R_P0(ndofs_0);
    // Minus(minus) projected_A_P0
    SimpleVector<double> projected_A_1(ndofs_1);
    matrix_A_1.multiply_with_vector(sln_1.v, projected_A_1.v, true);

    SimpleVector<double>* projected_A_P_0
      = (SimpleVector<double>*)cut_off_linear_part(projected_A_1.v, space_0, space_1);

    SimpleVector<double>* sln_1_projected = cut_off_linear_part(sln_1.v, space_0, space_1);
    matrix_A_0.multiply_with_vector(sln_1_projected->v, R_P0.v, true);

    SimpleVector<double> projected_f_P1(ndofs_1);
    projected_f_P1.set_vector(&f_P1);

    R_P0.change_sign();
    f_P0.add_vector(&R_P0);
    f_P0.add_vector(projected_A_P_0);
    f_P0.add_vector(cut_off_linear_part(projected_f_P1.v, space_0, space_1)->change_sign());
    f_P0.change_sign();

    num_coarse++;

    for(int iteration = 1; iteration <= 1000; iteration++)
    {
      SimpleVector<double>* rhs;
      if(iteration == 1)
      {
        if(polynomialDegree > 1)
        {
          rhs = (SimpleVector<double>*)cut_off_linear_part(residual_1, space_0, space_1)->add_vector(cut_off_linear_part(projected_f_P1.v, space_0, space_1)->change_sign());
          memcpy(vector_A_0.v, rhs->v, ndofs_0 * sizeof(double));
        }
        else
        {
          matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);
          vector_A_0.change_sign()->add_vector(&vector_b_0);
        }
      }
      else
      {
        // A(u_K) - done after the first step.
        matrix_A_0.multiply_with_vector(sln_0.v, vector_A_0.v, true);

        if(polynomialDegree > 1)
          vector_A_0.change_sign()->add_vector(&f_P0)->add_vector(&vector_b_0);
        else
          vector_A_0.change_sign()->add_vector(&vector_b_0);
      }

      solver_0.solve();
      sln_0.add_vector(solver_0.get_sln_vector());

      if(show_intermediate)
      {
        // Make solution
        Solution<double>::vector_to_solution(&sln_0, space_0, previous_sln);
        solution_view->show(previous_sln);
      }

      // Residual check.
      residual_condition(&matrix_A_0, &vector_b_0, sln_0.v, residual_0, logger, iteration, false);

      if(get_l2_norm(solver_0.get_sln_vector(), space_0->get_num_dofs()) < 1e-6)
        break;
    }
#pragma endregion

#pragma region 1 - intermediate level
    // Store the previous solution.
    sln_1.set_vector(merge_slns(sln_0.v, space_0, sln_1.v, space_1, space_1, false));
    Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);

    for(int iteration = 1; iteration <= steps; iteration++)
    {
      // Solve for increment.
      matrix_A_1.multiply_with_vector(sln_1.v, vector_A_1.v, true);

      if(polynomialDegree > 1)
        vector_A_1.change_sign()->add_vector(&f_P1)->add_vector(&vector_b_1);
      else
        vector_A_1.change_sign()->add_vector(&vector_b_1);

      solver_1.solve();
      sln_1.add_vector(solver_1.get_sln_vector());

      // Make solution
      Solution<double>::vector_to_solution(&sln_1, space_1, previous_sln);
      if(show_intermediate)
        solution_view->show(previous_sln);

      // Residual check.
      residual_condition(&matrix_A_1, &vector_b_1, sln_1.v, residual_1, logger, iteration, false);
    }
#pragma endregion

#pragma region 0 - highest level

    if(polynomialDegree > 1)
    {
      // Store the previous solution.
      sln_2.set_vector(merge_slns(sln_1.v, space_1, sln_2.v, space_2, space_2, false));
      Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);

      for(int iteration = 1; iteration <= steps; iteration++)
      {
        // Solve for increment.
        matrix_A_2.multiply_with_vector(sln_2.v, vector_A_2.v, true);
        vector_A_2.change_sign()->add_vector(&vector_b_2);
        solver_2.solve();
        sln_2.add_vector(solver_2.get_sln_vector());

        // Make solution
        Solution<double>::vector_to_solution(&sln_2, space_2, previous_sln);
        solution_view->show(previous_sln);

        // Residual check.
        residual_condition(&matrix_A_2, &vector_b_2, sln_2.v, residual_2, logger, iteration, true);
      }
    }
#pragma endregion

    // Error & exact solution display.
    solution_view->show(previous_sln);

    if(error_condition(calc_l2_error(solvedExample, mesh, previous_sln, es, logger)))
      break;
  }
  logger.info("V-cycles: %i", v_cycles);
  logger.info("Coarse systems solved: %i", num_coarse);
  logger.info("1 systems solved: %i", num_1);
  logger.info("2 systems solved: %i", num_2);
}
