#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;
// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;
int REFINEMENT_COUNT = 0;
int REFINEMENT_COUNT_THRESHOLD = 5;

// Shock capturing.
bool SHOCK_CAPTURING = true;
// Initial polynomial degree.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Initial time step.
double time_step_length = 1E-4;
double TIME_INTERVAL_LENGTH = .2;
// CFL value.
double CFL_NUMBER = 0.3; 
// Kappa.
const double KAPPA = 1.4;

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

class CustomErrorCalculator : public ErrorCalculator<double>
{
public:
  CustomErrorCalculator(CalculatedErrorType errorType, int component_count) : ErrorCalculator<double>(errorType)
  {
    for(int i = 0; i < component_count; i++)
    {
      this->add_error_form(new CustomNormFormVol(i, i));
    }
  }

  class CustomNormFormVol : public NormFormVol<double>
  {
  public:
    CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
    {
    }

    double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
    {
      double result = 0.;
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->val[i];
      return result * e->area;
    }
  };
};

// Stopping criterion for adaptivity.
bool adaptivityErrorStop(int iteration, double error, int ref_ndof)
{
  if(ref_ndof < 6e3)
    return false;
  if(ref_ndof > 10e3)
    return true;

  return(error < 1e-2);
}

int main(int argc, char* argv[])
{
HermesCommonApi.set_integral_param_value(numThreads, 1);
  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  Hermes::Mixins::Loggable logger(true);
  logger.set_logFile_name("computation.log");

#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Perform initial mesh refinements.
  mesh->refine_in_area("Pre", 1, 2, true);
  mesh->refine_in_area("Pre", 1, 2, true);
  mesh->refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);

  int ndof = Space<double>::get_num_dofs(spaces);
  logger.info("Ndof: %d", ndof);
#pragma endregion

#pragma region 1.1 Initialize solutions.
  // Set initial conditions.
  MeshFunctionSharedPtr<double> exact_rho(new CustomInitialCondition(mesh, 0, KAPPA));
  MeshFunctionSharedPtr<double> exact_rho_v_x(new CustomInitialCondition(mesh, 1, KAPPA));
  MeshFunctionSharedPtr<double> exact_rho_v_y(new CustomInitialCondition(mesh, 2, KAPPA));
  MeshFunctionSharedPtr<double> exact_e(new CustomInitialCondition(mesh, 3, KAPPA));
  MeshFunctionSharedPtr<double> prev_rho(new CustomInitialCondition(mesh, 0, KAPPA));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new CustomInitialCondition(mesh, 1, KAPPA));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new CustomInitialCondition(mesh, 2, KAPPA));
  MeshFunctionSharedPtr<double> prev_e(new CustomInitialCondition(mesh, 3, KAPPA));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  MeshFunctionSharedPtr<double> sln_rho(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_e(new Solution<double>);
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);
  MeshFunctionSharedPtr<double> rsln_rho(new CustomInitialCondition(mesh, 0, KAPPA));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new CustomInitialCondition(mesh, 1, KAPPA));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new CustomInitialCondition(mesh, 2, KAPPA));
  MeshFunctionSharedPtr<double> rsln_e(new CustomInitialCondition(mesh, 3, KAPPA));
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);
#pragma endregion

#pragma region 2 Initialize weakform.
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> prescribed_markers("BottomPost", "TopPre", "TopPost", "Left");
  solid_wall_markers.push_back("BottomRefl");

  EulerEquationsWeakFormExplicitDoubleReflection wf(KAPPA, solid_wall_markers, prescribed_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, exact_rho, exact_rho_v_x, exact_rho_v_y, exact_e, (P_INIT == 0));
#pragma endregion

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(rslns));

  ScalarView density_view("Density", new WinGeom(0, 0, 600, 300));
  ScalarView pressure_view("Pressure", new WinGeom(0, 330, 600, 300));
  ScalarView velocity_view("Velocity", new WinGeom(650, 0, 600, 300));
  ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

#pragma region Solver setup
  LinearSolver<double> solver;
  solver.set_weak_formulation(&wf);
  solver.set_verbose_output(false);
#pragma endregion

#pragma region Adaptivity setup
  // Initialize refinement selector.
  HOnlySelector<double> selector;
  Hermes::vector<RefinementSelectors::Selector<double> *> selectors(&selector, &selector, &selector, &selector);

  // Error calculation.
  CustomErrorCalculator errorCalculator(RelativeErrorToGlobalNorm, 4);
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionCumulative<double> stoppingCriterion(.75);
  Adapt<double> adaptivity(spaces, &errorCalculator, &stoppingCriterion);
#pragma endregion

#pragma region 3.1 Set up reference mesh and spaces.
  Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
  MeshSharedPtr ref_mesh;
  SpaceSharedPtr<double> ref_space_rho;
  SpaceSharedPtr<double> ref_space_rho_v_x;
  SpaceSharedPtr<double> ref_space_rho_v_y;
  SpaceSharedPtr<double> ref_space_e;
#pragma endregion

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    logger.info("Time step %d, time %3.5f, time step %3.5f.", iteration, t, time_step_length);

#pragma region 4.1. Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0) 
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      Space<double>::unrefine_all_mesh_elements(spaces);
      if(REFINEMENT_COUNT > REFINEMENT_COUNT_THRESHOLD)
        Space<double>::unrefine_all_mesh_elements(spaces);
      REFINEMENT_COUNT = 0;
      Space<double>::assign_dofs(spaces);
      solver.free_cache();
    }
#pragma endregion

#pragma region 4.2. Adaptivity loop.
    int as = 1;
    do
    {
#pragma region 7.1 Create reference mesh and spaces.
      ref_mesh = refMeshCreatorFlow.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, 0);
      ref_space_rho = refSpaceCreatorRho.create_ref_space();
      ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      ref_space_e = refSpaceCreatorE.create_ref_space();
      Hermes::vector<SpaceSharedPtr<double>  > ref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);
      solver.set_spaces(ref_spaces);
#pragma endregion

      int ref_ndof = Space<double>::get_num_dofs(ref_spaces);
      logger.info("\tAdaptivity step %d, NDOFs: %i.", as, ref_ndof);

#pragma region 7.2 Solve
      logger.info("\t\tSolving.");
      ((CustomInitialCondition*)exact_rho.get())->time = t;
      ((CustomInitialCondition*)exact_rho_v_x.get())->time = t;
      ((CustomInitialCondition*)exact_rho_v_y.get())->time = t;
      ((CustomInitialCondition*)exact_e.get())->time = t;

      // Set the current time step.
      wf.set_current_time_step(time_step_length);
      solver.solve();
#pragma endregion

#pragma region 7.3 Get the solution with optional shock capturing.
      if(!SHOCK_CAPTURING)
      {
        logger.info("\t\tObtaining solution.");
        Solution<double>::vector_to_solutions(solver.get_sln_vector(), ref_spaces, rslns);
      }
      else
      {
        logger.info("\t\tObtaining limited solution.");
        PostProcessing::VertexBasedLimiter limiter(ref_spaces, solver.get_sln_vector(), 1);
        limiter.get_solutions(rslns);
        limitVelocityAndEnergy(ref_spaces, limiter, rslns);
      }

      // Calculate time step according to CFL condition.
      CFL.calculate(rslns, (ref_spaces)[0]->get_mesh(), time_step_length);
#pragma endregion

#pragma region 7.3.1 Visualization
      if((iteration % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION)
        {        
          //pressure->reinit();
          velocity->reinit();
          density_view.show(rslns[0]);
          //pressure_view.show(pressure);
          //velocity_view.show(velocity);
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION)
        {
          pressure->reinit();
          velocity->reinit();
          Linearizer lin;
          char filename[40];
          sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          lin.save_solution_vtk(pressure, filename, "Pressure", false);
          sprintf(filename, "Velocity-%i.vtk", iteration - 1);
          lin.save_solution_vtk(velocity, filename, "Velocity", false);
          sprintf(filename, "Rho-%i.vtk", iteration - 1);
          lin.save_solution_vtk(prev_rho, filename, "Rho", false);
        }
      }
#pragma endregion


#pragma region 7.4 Project to coarse mesh -> error estimation -> space adaptivity
      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("\t\tProjecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(spaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(slns, rslns);
      double density_error = errorCalculator.get_error_squared(0) * 100;

      if((iteration % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
      {
        if(HERMES_VISUALIZATION)
          eview.show(errorCalculator.get_errorMeshFunction(0));
      }

      // Report results.
      logger.info("\t\tDensity error: %g%%.", density_error);

      // If err_est too large, adapt the mesh.
      if (adaptivityErrorStop(iteration, density_error, ref_ndof))
        break;
      else
      {
        Hermes::Mixins::Loggable::Static::info("\t\tAdapting coarse mesh.");
        adaptivity.adapt(selectors);
        REFINEMENT_COUNT++;
        as++;
      }
#pragma endregion
    }
    while (true);
#pragma endregion

    // Copy the solutions into the previous time level ones.
    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);

    iteration++;
  }
#pragma endregion

  // Done.
  return 0;
}
