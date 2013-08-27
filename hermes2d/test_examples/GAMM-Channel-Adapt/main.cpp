#include "hermes2d.h"
#include "../euler_util.h"

double CFL_NUMBER = 1.0;

#pragma region Visulization + Utilities
double LAST_UNREF = 0.;
// Initial time step.
double time_step_length = 1e-6;
double TIME_INTERVAL_LENGTH = 20.;

// Kappa.
const double KAPPA = 1.4;

// Set up CFL calculation class.
CFLCalculation CFL(CFL_NUMBER, KAPPA);

// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.5;         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.25;       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        

// Mesh filename.
const std::string MESH_FILENAME = "GAMM-channel.mesh";

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";
#pragma endregion

#pragma region Limiting
// Limiting
bool SHOCK_CAPTURING = false;
const EulerLimiterType limiter_type = VertexBased;
bool limit_velocities = false;
#pragma endregion

// Initial polynomial degree.
const int P_INIT = 1;
class CustomErrorCalculator : public ErrorCalculator<double>
{
public:
  CustomErrorCalculator(CalculatedErrorType errorType, int component_count) : ErrorCalculator<double>(errorType)
  {
    for(int i = 0; i < component_count; i++)
    {
      //this->add_error_form(new CustomNormFormVol(i, i));
      this->add_error_form(new CustomNormFormDG(i, i));
    }
  }

  class CustomNormFormVol : public NormFormVol<double>
  {
  public:
    CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
    {
      //this->functionType = CoarseSolutions;
    }

    double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
    {
      double result_values = 0., result_derivatives = 0.;
      for (int point_i = 0; point_i < n; point_i++)
      {
        result_values += wt[point_i] * (u->val[point_i] * v->val[point_i]);
        result_derivatives += wt[point_i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }
      return result_values;// + result_derivatives;
    }
  };

  class CustomNormFormDG : public NormFormDG<double>
  {
  public:
    CustomNormFormDG(int i, int j) : NormFormDG<double>(i, j)
    {
      this->functionType = CoarseSolutions;
    }

    double value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
    {
      double result = double(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]);
      return result;
    }
  };
};

bool adaptivityErrorStop(int iteration, double time, double error, int ref_ndof)
{
  if(iteration == 0)
    return true;
  return false;
  return error < .001;
}

// Number of initial uniform mesh refinements.
int INIT_REF_NUM = 2;

// Every UNREF_FREQth time step the mesh is unrefined.
double UNREF_FREQ = 1e-2;

// Error calculation.
CustomErrorCalculator errorCalculator(AbsoluteError, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionCumulative<double> stoppingCriterion(.6);
// Initialize refinement selector.
HOnlySelector<double> selector;
int order_increase = 0;
int max_p = 0;

int main(int argc, char* argv[])
{
  Hermes::Mixins::Loggable logger(true);
  logger.set_logFile_name("computation.log");

#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh.
  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> inlet_markers;
  Hermes::vector<std::string> outlet_markers;
  MeshSharedPtr mesh(new Mesh);

  MeshReaderH2D mloader;
  mloader.load(MESH_FILENAME.c_str(), mesh);
  solid_wall_markers.push_back(BDY_SOLID_WALL_BOTTOM);
  solid_wall_markers.push_back(BDY_SOLID_WALL_TOP);
  inlet_markers.push_back(BDY_INLET);
  outlet_markers.push_back(BDY_OUTLET);

  // Perform initial mesh refinements.
  mesh->refine_element_id(1, 2);
  mesh->refine_all_elements(1);
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements(0, true);

  mesh->refine_towards_vertex(1, 2);
  mesh->refine_towards_vertex(2, 2);
    mesh->refine_all_elements(0, true);
  
  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);

  SpaceSharedPtr<double> const_space_rho(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_rho_v_x(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_rho_v_y(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  SpaceSharedPtr<double> const_space_e(new L2Space<double>(mesh, 0, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > const_spaces(const_space_rho, const_space_rho_v_x, const_space_rho_v_y, const_space_e);

  int ndof = Space<double>::get_num_dofs(spaces);
  logger.info("Ndof: %d", ndof);
#pragma endregion

#pragma region 2. Prev slns
  // Set initial conditions.
  MeshFunctionSharedPtr<double> sln_rho(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double>);
  MeshFunctionSharedPtr<double> sln_e(new Solution<double>);
  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);

  MeshFunctionSharedPtr<double> rsln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> rsln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

  // For updating of mean values.
  MeshFunctionSharedPtr<double> updated_prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> updated_prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> updated_prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> updated_prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
  Hermes::vector<MeshFunctionSharedPtr<double> > updated_prev_slns(updated_prev_rho, updated_prev_rho_v_x, updated_prev_rho_v_y, updated_prev_e);
#pragma endregion

  EulerEquationsWeakFormSemiImplicit wf_implicit(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,
    solid_wall_markers, inlet_markers, outlet_markers,
    prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, true);

  EulerEquationsWeakFormExplicit wf_explicit(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,
    solid_wall_markers, inlet_markers, outlet_markers,
    prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, updated_prev_rho, updated_prev_rho_v_x, updated_prev_rho_v_y, updated_prev_e, (P_INIT == 0));

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(rslns));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(650, 0, 600, 300));
  ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

#pragma region Solver setup
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  LinearSolver<double> solver_implicit(&wf_implicit, const_spaces);
  ((IterSolver<double>*)solver_implicit.get_linear_solver())->set_tolerance(1e-6, IterSolver<double>::AbsoluteTolerance);
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_UMFPACK);
  LinearSolver<double> solver_explicit(&wf_explicit, spaces);
  solver_explicit.set_jacobian_constant();
  wf_explicit.set_current_time_step(time_step_length);
  wf_implicit.set_current_time_step(time_step_length);
  DiscreteProblemDGAssembler<double>::dg_order = 10;
#pragma endregion

#pragma region Adaptivity setup
  Hermes::vector<RefinementSelectors::Selector<double> *> selectors(&selector, &selector, &selector, &selector);
  Adapt<double> adaptivity(space_rho, &errorCalculator, &stoppingCriterion);
#pragma endregion

#pragma region 3.1 Set up reference mesh and spaces.
  Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
  MeshSharedPtr ref_mesh;
  SpaceSharedPtr<double> ref_space_rho;
  SpaceSharedPtr<double> ref_space_rho_v_x;
  SpaceSharedPtr<double> ref_space_rho_v_y;
  SpaceSharedPtr<double> ref_space_e;
  Hermes::vector<SpaceSharedPtr<double>  > ref_spaces;

  SpaceSharedPtr<double> const_ref_space_rho;
  SpaceSharedPtr<double> const_ref_space_rho_v_x;
  SpaceSharedPtr<double> const_ref_space_rho_v_y;
  SpaceSharedPtr<double> const_ref_space_e;
  Hermes::vector<SpaceSharedPtr<double> > const_ref_spaces;
#pragma endregion

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    logger.info("Time step %d, time %3.5f, time step %3.5f.", iteration, t, time_step_length);

#pragma region 4.1. Periodic global derefinements.
    if (iteration > 1 && t > LAST_UNREF + UNREF_FREQ)
    {
      LAST_UNREF = t;
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      space_rho->unrefine_all_mesh_elements(true);
      space_rho->adjust_element_order(-1, P_INIT);
      space_rho_v_x->adjust_element_order(-1, P_INIT);
      space_rho_v_y->adjust_element_order(-1, P_INIT);
      space_e->adjust_element_order(-1, P_INIT);

      const_space_rho->set_uniform_order(0);
      const_space_rho_v_x->set_uniform_order(0);
      const_space_rho_v_y->set_uniform_order(0);
      const_space_e->set_uniform_order(0);

      Space<double>::assign_dofs(spaces);
    }
#pragma endregion

#pragma region 4.2. Adaptivity loop.
    int as = 1;
    do
    {
#pragma region 7.1 Create reference mesh and spaces.
      ref_mesh = refMeshCreatorFlow.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, order_increase);
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, order_increase);

      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRho(const_space_rho, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRhoVx(const_space_rho_v_x, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorRhoVy(const_space_rho_v_y, ref_mesh, 0);
      Space<double>::ReferenceSpaceCreator const_refSpaceCreatorE(const_space_e, ref_mesh, 0);

      ref_space_rho = refSpaceCreatorRho.create_ref_space();
      ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      ref_space_e = refSpaceCreatorE.create_ref_space();
      ref_spaces = Hermes::vector<SpaceSharedPtr<double> >(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);
      solver_explicit.set_spaces(ref_spaces);

      const_ref_space_rho = const_refSpaceCreatorRho.create_ref_space();
      const_ref_space_rho_v_x = const_refSpaceCreatorRhoVx.create_ref_space();
      const_ref_space_rho_v_y = const_refSpaceCreatorRhoVy.create_ref_space();
      const_ref_space_e = const_refSpaceCreatorE.create_ref_space();
      const_ref_spaces = Hermes::vector<SpaceSharedPtr<double> >(const_ref_space_rho, const_ref_space_rho_v_x, const_ref_space_rho_v_y, const_ref_space_e);
      solver_implicit.set_spaces(const_ref_spaces);
#pragma endregion

      int ref_ndof = Space<double>::get_num_dofs(ref_spaces);
      logger.info("\tAdaptivity step %d, NDOFs: %i.", as, ref_ndof);

      if(iteration)
      {
        // Solve.
        // 0 - store the sln vector.
        double* previous_sln_vector = new double[ref_ndof];
        OGProjection<double>::project_global(ref_spaces, prev_slns, previous_sln_vector);

        // 1 - solve implicit.
        solver_implicit.solve();
        double* mean_values = new double[ref_ndof];
        Solution<double>::vector_to_solutions(solver_implicit.get_sln_vector(), const_ref_spaces, prev_slns);
        OGProjection<double>::project_global(ref_spaces, prev_slns, mean_values);

        // 2 - Update the mean values.
        Space<double>::assign_dofs(ref_spaces);
        for(int component = 0; component < 4; component++)
        {
          Element* e;
          for_all_active_elements(e, ref_mesh)
          {
            AsmList<double> al_fine;
            ref_spaces[component]->get_element_assembly_list(e, &al_fine);
            for(unsigned int shape_i = 0; shape_i < al_fine.cnt; shape_i++)
            {
              int order = ref_spaces[component]->get_shapeset()->get_order(al_fine.idx[shape_i], e->get_mode());
              if(order == 0)
              {
                int dof = al_fine.dof[shape_i];
                previous_sln_vector[dof] = mean_values[dof];
              }
            }
          }
        }

        // Solve explicit.
        Solution<double>::vector_to_solutions(previous_sln_vector, ref_spaces, prev_slns);

        // Clean up.
        delete [] previous_sln_vector;
        delete [] mean_values;
      }

      solver_explicit.solve();

#pragma region *. Get the solution with optional shock capturing.
      if(!SHOCK_CAPTURING)
        Solution<double>::vector_to_solutions(solver_explicit.get_sln_vector(), ref_spaces, rslns);
      else
      {
        PostProcessing::Limiter<double>* limiter = create_limiter(limiter_type, ref_spaces, solver_explicit.get_sln_vector(), max_p);
        limiter->get_solutions(rslns);
        if(limit_velocities)
          limitVelocityAndEnergy(ref_spaces, limiter, rslns);
        delete limiter;
      }
#pragma endregion

#pragma region 7.4 Project to coarse mesh -> error estimation -> space adaptivity
      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("\t\tProjecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(spaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(sln_rho, rsln_rho);
      double density_error = errorCalculator.get_error_squared(0) * 100;
      if(iteration)
      {
        std::cout << errorCalculator.get_element_error_squared(0, 38) << std::endl;
        std::cout << errorCalculator.get_element_error_squared(0, 21) << std::endl;
      }
      if(HERMES_VISUALIZATION)
        eview.show(errorCalculator.get_errorMeshFunction(0));

      // Report results.
      logger.info("\t\tDensity error: %g%%.", density_error);

#pragma region 7.3.1 Visualization
      if((iteration % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION)
        {        
          pressure->reinit();
          Mach_number->reinit();
          pressure_view.show(slns[0]);
          Mach_number_view.show(Mach_number);
          View::wait_for_keypress();
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION)
        {
          pressure->reinit();
          Mach_number->reinit();
          Linearizer lin;
          char filename[40];
          sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          lin.save_solution_vtk(pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(Mach_number, filename, "Velocity", false);
          Orderizer ord;
          sprintf(filename, "Mesh-%i.vtk", iteration - 1);
          ord.save_mesh_vtk(ref_space_rho, filename);
        }
      }
#pragma endregion

      // If err_est too large, adapt the mesh.
      if (adaptivityErrorStop(iteration, t, density_error, ref_ndof))
        break;
      else
      {
        Hermes::Mixins::Loggable::Static::info("\t\tAdapting coarse mesh.");
        adaptivity.adapt(&selector);
        space_rho_v_x->copy(space_rho, mesh);
        space_rho_v_y->copy(space_rho, mesh);
        space_e->copy(space_rho, mesh);
        Space<double>::assign_dofs(spaces);

        const_space_rho->copy(space_rho, mesh);
        const_space_rho_v_x->copy(space_rho, mesh);
        const_space_rho_v_y->copy(space_rho, mesh);
        const_space_e->copy(space_rho, mesh);

        const_space_rho->set_uniform_order(0);
        const_space_rho_v_x->set_uniform_order(0);
        const_space_rho_v_y->set_uniform_order(0);
        const_space_e->set_uniform_order(0);

        Space<double>::assign_dofs(const_spaces);
        as++;
      }
#pragma endregion
    }
    while (true);

    // Calculate time step according to CFL condition.
    CFL.calculate(solver_explicit.get_sln_vector(), ref_spaces, time_step_length);

    wf_explicit.set_current_time_step(time_step_length);
    wf_implicit.set_current_time_step(time_step_length);

    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);

    updated_prev_rho->copy(rsln_rho);
    updated_prev_rho_v_x->copy(rsln_rho_v_x);
    updated_prev_rho_v_y->copy(rsln_rho_v_y);
    updated_prev_e->copy(rsln_e);

    iteration++;
  }
#pragma endregion

  // Done.
  View::wait();
  return 0;
}
