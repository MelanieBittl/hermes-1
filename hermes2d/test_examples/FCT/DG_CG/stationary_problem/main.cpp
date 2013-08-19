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
bool SD = false;

bool serendipity = false;

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


//perturbed mesh

/*Node* vn; Element* e; double h; 
for_all_active_elements(e, mesh) 
{
	h = Hermes::sqrt(0.5)*e->get_diameter(); // Laenge der kleinste Kante
	break;
}
for_all_vertex_nodes(vn, mesh)
{
	double delta_x = (rand()%50)/100.; 
	double delta_y = (rand()%50)/100.; 
		if(!vn->bnd) 
		{		vn->x +=h*(delta_x-0.5) ;
				vn->y +=h*(delta_y-0.5) ;
		}
}
*/


  // Create an space with default shapeset.
//CustomDirichletCondition bc_essential(Hermes::vector<std::string>("inlet1","inlet2"));
 // EssentialBCs<double>  bcs(&bc_essential);
  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT, serendipity));	
 //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));


/*
	Orderizer ord;
	ord.save_mesh_vtk(space, "perturbed_mesh_50_ref4.vtk");

MeshView mview("Meshview", new WinGeom(450, 0, 440, 350));
mview.show(mesh);
View::wait(HERMES_WAIT_KEYPRESS);


BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
bview.show(space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));
 
  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView s2view("solution_2", new WinGeom(0, 500, 500, 400));

	ScalarView fview("filter", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	//lview.show(u_prev_time);
//fview.set_min_max_range(0,0.001);
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

	CustomWeakForm wf_surf(u_prev_time,mesh,all, DG,SD);
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

/*
FILE * matFile;
matFile = fopen ("test.m","w");
if(dg_surface_matrix->dump(matFile, "mat_t")) printf("print matrix\n");
fclose (matFile);  

/*
CustomWeakForm wf(u_prev_time,mesh,true, false);
 // Initialize linear solver.
 Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();
	for(int i=0; i<ndof; i++) coeff_vec_2[i] = sln_vector[i];

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, u_new);
		s2view.show(u_new);
  }catch(std::exception& e)
  {
    std::cout << e.what();
  }
*/

calc_error_total(u_new, u_prev_time,space);






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

