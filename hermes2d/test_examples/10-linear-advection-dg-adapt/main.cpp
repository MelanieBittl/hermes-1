#include "definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//  It is intended to show how evalutation of surface matrix forms that take basis functions defined
//  on different elements work. It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:    Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on[0,0.5] x {0}, and g = 0 anywhere else.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF = 3;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 1;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.5;

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_H_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e11;

int main(int argc, char* args[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++)
    mesh->refine_all_elements();

  // Create an L2 space.
  SpaceSharedPtr<double> coarse_space(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> fine_space(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST);
  selector.set_error_weights(1.,1.,1.);

  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> ref_coarse_sln(new Solution<double>);
  MeshFunctionSharedPtr<double> ref_fine_sln(new Solution<double>);

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_bottom_left", mesh);
  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver;
  linear_solver.set_weak_formulation(&wf);

  adaptivity.set_space(fine_space);

  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh
    // and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator ref_coarse_space_creator(coarse_space, ref_mesh, 0);
    SpaceSharedPtr<double> ref_coarse_space = ref_coarse_space_creator.create_ref_space();

    Space<double>::ReferenceSpaceCreator ref_fine_space_creator(fine_space, ref_mesh, 0);
    SpaceSharedPtr<double> ref_fine_space = ref_fine_space_creator.create_ref_space();

    // Solve the linear system. If successful, obtain the solution.
#pragma region old_implementation
    /*
    try
    {
      linear_solver.set_space(ref_coarse_space);
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_coarse_space, ref_coarse_sln);

      view1.show(ref_coarse_sln);
      view1.set_title("Coarse solution");
      view1.wait_for_keypress();

      linear_solver.set_space(ref_fine_space);
      linear_solver.solve();

      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_fine_space, ref_fine_sln);
      view1.show(ref_fine_sln);
      view1.set_title("Unlimited fine solution");
      view1.wait_for_keypress();

      FluxLimiter* flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, linear_solver.get_sln_vector(), ref_fine_space, true);
      flux_limiter->limit_according_to_detector();
      flux_limiter->get_limited_solution(ref_fine_sln);

      view1.show(ref_fine_sln);
      view1.set_title("Limited fine solution - UNseparated coarse / fine component");
      view1.wait_for_keypress();

      Element* e;
      for_all_active_elements(e, ref_fine_sln->get_mesh())
      {
        AsmList<double> al;
        ref_fine_space->get_element_assembly_list(e, &al);
        for(unsigned int shape_i = 0; shape_i < al.get_cnt(); shape_i++)
          if(ref_fine_space->get_shapeset()->get_order(al.get_idx()[shape_i], e->get_mode()) == 0)
            linear_solver.get_sln_vector()[al.get_dof()[shape_i]] = 0.0;
      }

      double* coarse_sln_vector = new double[ref_fine_space->get_num_dofs()];
      OGProjection<double>::project_global(ref_fine_space, ref_coarse_sln, coarse_sln_vector);
      for(int i = 0; i < ref_fine_space->get_num_dofs(); i++)
        coarse_sln_vector[i] += linear_solver.get_sln_vector()[i];

      Solution<double>::vector_to_solution(coarse_sln_vector, ref_fine_space, ref_fine_sln);

      view1.set_title("Limited fine solution - separated coarse / fine component");
      view1.show(ref_fine_sln);
      OGProjection<double>::project_global(fine_space, ref_fine_sln, sln, HERMES_L2_NORM);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
    */
#pragma endregion

#pragma region new_implementation
    try
    {
      linear_solver.set_space(ref_fine_space);
      linear_solver.solve();

      PostProcessing::VertexBasedLimiter limiter(ref_fine_space, linear_solver.get_sln_vector());
      ref_fine_sln = limiter.get_solution();

      // Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_fine_space, ref_fine_sln);
      view1.show(ref_fine_sln);
      OGProjection<double>::project_global(fine_space, ref_fine_sln, sln, HERMES_L2_NORM);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
#pragma endregion

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_fine_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    std::cout << "Error: " << err_est_rel << "%." << std::endl;

    // If err_est_rel too large, adapt the mesh->
    if(err_est_rel < ERR_STOP)
      done = true;
    else
      done = adaptivity.adapt(&selector);
    as++;
  }
  while (done == false);

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}
