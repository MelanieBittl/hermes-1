#define HERMES_REPORT_ALL
#include "definitions.h"
#include "hermes2d.h"
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

const int INIT_REF_NUM =8;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

bool all = true;
bool DG = false;
bool SD = true;


#include "error_estimates.cpp"

int main(int argc, char* argv[])
{ 
  // Load the mesh..
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements(); 
/*MeshView m1view("Meshview1", new WinGeom(450, 450, 440, 350));
m1view.show(mesh);*/

  // Create an space with default shapeset.
  //SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
 //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));


  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));
 
  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	//lview.show(u_prev_time);
	
	int ndof = space->get_num_dofs();

	double *coeff_vec = new double[ndof];		
	double *coeff_vec_2 = new double[ndof];
	double *coeff_vec_3 = new double[ndof];
	double *coeff_vec_4 = new double[ndof];
	double* vec_new;
  memset(coeff_vec_2, 0, ndof*sizeof(double));
 memset(coeff_vec_3, 0, ndof*sizeof(double));
 memset(coeff_vec_4, 0, ndof*sizeof(double));


	///////////////////////////------------false, false, false (only CG)-----------
//																			all, 		DG, 		SD ---------------------------------------------------------------------

	CustomWeakForm wf_surf(u_prev_time,mesh,all, DG,SD);
	CSCMatrix<double>* dg_surface_matrix = new CSCMatrix<double> ; 
	SimpleVector<double> * surf_rhs = new SimpleVector<double> (ndof); 
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,space);	
	//dp_surf->set_linear(true,false);
	dp_surf->assemble(dg_surface_matrix,surf_rhs);
	UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>( dg_surface_matrix, surf_rhs);    
	try
	{
	solver->solve();
	}
	catch(Hermes::Exceptions::Exception e)
	{
		e.print_msg();
	}	
	vec_new = solver->get_sln_vector();
	Solution<double>::vector_to_solution(vec_new, space, u_new);
	for(int i=0; i<ndof; i++) coeff_vec_2[i] = vec_new[i];

		//sview.show(u_new);

calc_error_total(u_new, u_prev_time,space);

	
  // Wait for the view to be closed.
  View::wait();
  return 0;
}

