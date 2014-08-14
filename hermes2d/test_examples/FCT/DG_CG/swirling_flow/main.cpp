#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM = 3;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.
                     
const double time_step = 2e-3;                           // Time step.                      
const double T_FINAL = 1.5;  

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)


const bool all = true;
const bool DG = true;
const bool SD = false;

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";


#include "error_estimations.cpp"
#include "mass_lumping.cpp"


int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
mloader.load("unit.mesh", basemesh);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);


  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
 //SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));	
  
  int ndof = space->get_num_dofs();


  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition (mesh));
  MeshFunctionSharedPtr<double>  u_proj(new Solution<double>);

// Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;
//lin.save_solution_vtk(u_prev_time, "init.vtk", "solution", mode_3D);
 
  // Initialize views.
	ScalarView sview("Loesung", new WinGeom(500, 500, 500, 400));
	ScalarView lview("init Loesung", new WinGeom(500, 0, 500, 400));
	//lview.show(u_prev_time);

  
	char title[100];
    int ts = 1;
    double current_time =0.;
	 
CustomWeakFormMassmatrix  massmatrix;
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;  
			DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, space);
			dp_mass->assemble(mass_matrix); 	
			
			//CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);
		double* coeff_vec = new double[ndof]; 
		OGProjection<double> ogProjection;	
		ogProjection.project_global(space, u_prev_time, coeff_vec,  HERMES_L2_NORM);	
		
				//Lumped_Projection::project_lumped(space, u_prev_time, coeff_vec,lumped_matrix);
		
      Solution<double>::vector_to_solution(coeff_vec, space, u_new);
		Solution<double>::vector_to_solution(coeff_vec, space, u_proj);
		delete [] coeff_vec;


CustomWeakForm wf(u_prev_time, mesh,T_FINAL,time_step, theta, all, DG, SD, true);
		//CustomWeakForm wf(u_new, mesh,T_FINAL,time_step, theta, all, DG, SD, true);
 // Initialize linear solver.
 Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  do
  {
	  Hermes::Mixins::Loggable::Static::info("time=%f, ndof = %i ", current_time,ndof); 
	  wf.set_current_time(current_time);
					// Solve the linear problem.
					try
					{
						linear_solver.solve();

						// Get the solution vector.
						double* sln_vector = linear_solver.get_sln_vector();

						// Translate the solution vector into the previously initialized Solution.
						Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, u_new);
						//sview.show(u_new);
 						if(ts==1) wf.set_ext(u_new);
						//if(current_time >= 0.75) View::wait(HERMES_WAIT_KEYPRESS);
					}catch(std::exception& e)
					{
						std::cout << e.what();
					}
	 	current_time += time_step;
	 	ts++;
	  }
  while (current_time<T_FINAL);


calc_error_total(u_new, u_proj,u_prev_time, space, T_FINAL, current_time);


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

