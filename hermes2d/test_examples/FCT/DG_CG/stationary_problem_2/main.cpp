#define HERMES_REPORT_ALL
#include "definitions.h"
#include "hermes2d.h"
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.



MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";


bool all = true;
bool DG = true;


bool serendipity = true;

#include "error_estimates.cpp"


int main(int argc, char* argv[])
{ 
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements(); 
/*MeshView m1view("Meshview1", new WinGeom(450, 450, 440, 350));
m1view.show(mesh);*/

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));




//CustomDirichletCondition bc_essential(Hermes::vector<std::string>("inlet1","inlet2"));
//EssentialBCNonConst bc_nonconst("inlet",u_prev_time); 
//EssentialBCs<double>  bcs(&bc_nonconst);
  // Create an space with default shapeset.
  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
 //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh,P_INIT));



//SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	

//SpaceSharedPtr<double> space_1(new H1Space<double>(mesh,P_INIT));	

// MeshFunctionSharedPtr<double>  u_test(new Solution<double>);
//dynamic_cast<Solution<double>*>(u_test.get())->load("solution", space);

  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	lview.show(u_prev_time);
//sview.set_min_max_range(0,1);
//View::wait(HERMES_WAIT_KEYPRESS);
	
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

	CustomWeakForm wf_surf(u_prev_time,mesh,all, DG);
	UMFPackMatrix<double>* dg_surface_matrix = new UMFPackMatrix<double> ; 
	UMFPackVector<double> * surf_rhs = new UMFPackVector<double> (ndof); 
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,space);	
	dp_surf->set_linear(true,false);
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

		sview.show(u_new);




calc_error_total(u_new, u_prev_time,space);






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

