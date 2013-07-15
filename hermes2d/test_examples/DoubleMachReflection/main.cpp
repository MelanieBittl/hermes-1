#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// Discontinuous Galerkin method of higher order with no adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Solid walls, inlet, outlet.
//
// IC: Constant state identical to inlet.
//

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

bool SHOCK_CAPTURING = false;

// Initial polynomial degree.
const int P_INIT = 0;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 6;
// CFL value.
double CFL_NUMBER = 0.1;
// Initial time step.
double time_step_length = 5E-4;
double TIME_INTERVAL_LENGTH = .2;

// Kappa.
const double KAPPA = 1.4;

// Weak forms.
#include "forms_explicit.cpp"

// Initial condition.
#include "initial_condition.cpp"

int main(int argc, char* argv[])
{
#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
  {
    mesh->refine_in_area("Post", 1, 1);
    mesh->refine_in_area("Pre", 1, 0);
  }

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);

  int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("Ndof: %d", ndof);
#pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

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

  Hermes::vector<std::string> solid_wall_markers;
  Hermes::vector<std::string> prescribed_markers("BottomPost", "TopPre", "TopPost", "Left");
  solid_wall_markers.push_back("BottomRefl");

  EulerEquationsWeakFormExplicitDoubleReflection wf(KAPPA, solid_wall_markers, prescribed_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, exact_rho, exact_rho_v_x, exact_rho_v_y, exact_e, (P_INIT == 0));

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double>  velocity(new VelocityFilter(prev_slns));

  ScalarView density_view("Density", new WinGeom(0, 0, 600, 300));
  ScalarView pressure_view("Pressure", new WinGeom(0, 330, 600, 300));
  ScalarView velocity_view("Velocity", new WinGeom(650, 0, 600, 300));
  ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
  ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

  LinearSolver<double> solver(&wf, spaces);
  solver.set_jacobian_constant();

#pragma region 4. Time stepping loop.
  int iteration = 0;
  for(double t = 0.0; t <= TIME_INTERVAL_LENGTH + Hermes::Epsilon; t += time_step_length)
  {
    // Info.
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_current_time_step(time_step_length);

    // Solve.
    ((CustomInitialCondition*)exact_rho.get())->time = t;
    ((CustomInitialCondition*)exact_rho_v_x.get())->time = t;
    ((CustomInitialCondition*)exact_rho_v_y.get())->time = t;
    ((CustomInitialCondition*)exact_e.get())->time = t;
    solver.solve();

#pragma region *. Get the solution with optional shock capturing.
    if(!SHOCK_CAPTURING)
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, prev_slns);
    else
    {
      PostProcessing::VertexBasedLimiter limiter(spaces, solver.get_sln_vector(), 1);
      limiter.get_solutions(prev_slns);

      int running_dofs = 0;
      int ndof = spaces[0]->get_num_dofs();
      double* density_sln_vector = limiter.get_solution_vector();
      Element* e;
      AsmList<double> al_density;
      for(int component = 1; component < 4; component++)
      {
        if(spaces[component]->get_num_dofs() != ndof)
          throw Exceptions::Exception("Euler code is supposed to be executed on a single mesh.");

        double* conservative_vector = limiter.get_solution_vector() + component * ndof;
        double* real_vector = new double[ndof];
        memset(real_vector, 0, sizeof(double) * ndof);

        for_all_active_elements(e, spaces[0]->get_mesh())
        {
          spaces[0]->get_element_assembly_list(e, &al_density);

          real_vector[al_density.dof[0]] = conservative_vector[al_density.dof[0]] / density_sln_vector[al_density.dof[0]];
          real_vector[al_density.dof[1]] = (conservative_vector[al_density.dof[1]] - real_vector[al_density.dof[0]] * density_sln_vector[al_density.dof[1]]) / density_sln_vector[al_density.dof[0]];
          real_vector[al_density.dof[2]] = (conservative_vector[al_density.dof[2]] - real_vector[al_density.dof[0]] * density_sln_vector[al_density.dof[2]]) / density_sln_vector[al_density.dof[0]];
        }

        PostProcessing::VertexBasedLimiter real_component_limiter(spaces[0], real_vector, 1);
        real_component_limiter.get_solution();
        real_vector = real_component_limiter.get_solution_vector();

        for_all_active_elements(e, spaces[0]->get_mesh())
        {
          spaces[0]->get_element_assembly_list(e, &al_density);

          conservative_vector[al_density.dof[1]] = density_sln_vector[al_density.dof[0]] * real_vector[al_density.dof[1]]
          + density_sln_vector[al_density.dof[1]] * real_vector[al_density.dof[0]];

          conservative_vector[al_density.dof[2]] = density_sln_vector[al_density.dof[0]] * real_vector[al_density.dof[2]]
          + density_sln_vector[al_density.dof[2]] * real_vector[al_density.dof[0]];
        }

        Solution<double>::vector_to_solution(conservative_vector, spaces[0], prev_slns[component]);
      }
    }
#pragma endregion

    // Calculate time step according to CFL condition.
    //CFL.calculate(prev_slns, mesh, time_step_length);

#pragma region 4.1. Visualization
    if(((iteration - 1) % EVERY_NTH_STEP == 0) || (t > TIME_INTERVAL_LENGTH - (time_step_length + Hermes::Epsilon)))
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        Mach_number->reinit();
        pressure->reinit();
        velocity->reinit();
        density_view.show(prev_slns[0]);
        pressure_view.show(pressure);
        velocity_view.show(velocity);
        order_view.show(space_rho);
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
  }
#pragma endregion

  // Done.
  View::wait();
  return 0;
}
