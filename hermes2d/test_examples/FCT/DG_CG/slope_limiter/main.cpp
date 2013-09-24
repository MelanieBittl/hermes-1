#define HERMES_REPORT_ALL
#include "definitions.h"
#include "shapeset_taylor.h" 
 #include "solution_slopelimiter.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 


const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
                    

  

bool all = true;
bool DG = true;
bool SD = false;


bool serendipity = true;

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


int main(int argc, char* argv[])
{	


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);
  
  

 // Initialize solution of lower & higher order
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double>  u_taylor(new Solution<double>);

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));	 


  // Initialize views.
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	ScalarView pview("sln", new WinGeom(500, 0, 500, 400));


  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
Orderizer ord;
	char title[100];



  OGProjection<double> ogProjection;

	SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
		
//SpaceSharedPtr<double> space_l2(new L2Space<double>(mesh,P_INIT));	

	int ndof = space->get_num_dofs();


	double *coeff_vec = new double[ndof];		
	double *coeff_vec_2 = new double[ndof];
	double *coeff_vec_3 = new double[ndof];
	double *coeff_vec_4 = new double[ndof];
	double* vec_new;
  memset(coeff_vec_2, 0, ndof*sizeof(double));
 memset(coeff_vec_3, 0, ndof*sizeof(double));
 memset(coeff_vec_4, 0, ndof*sizeof(double));

	CustomWeakForm wf_surf(u_prev_time,mesh,all, DG, SD);
	CSCMatrix<double>* dg_surface_matrix = new CSCMatrix<double> ; 
	SimpleVector<double> * surf_rhs = new SimpleVector<double> (ndof); 
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,space);	
	dp_surf->set_linear(true,true);
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



 SpaceSharedPtr<double> space_taylor(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
//SpaceSharedPtr<double> space_taylor(new L2Space<double>(mesh, P_INIT));
	double *coeff_vec_taylor = new double[space_taylor->get_num_dofs()];
ogProjection.project_global(space_taylor,u_new, coeff_vec_taylor,  HERMES_L2_NORM); 
Solution<double>::vector_to_solution(coeff_vec_taylor, space_taylor, u_taylor);



/*	MeshFunctionSharedPtr<double> sln(new SlopeLimiterSolution(u_taylor, space_taylor));	
(dynamic_cast<SlopeLimiterSolution*>(sln.get()))->limit_solution_according_to_detector();
//		sprintf(title, "lim. Loesung");
			//sview.set_title(title);
			sview.show(sln);*/




    PostProcessing::VertexBasedLimiter limiter_3(space_taylor, coeff_vec_taylor, 2.);
   u_new->copy(limiter_3.get_solution());


		sprintf(title, "lim. Loesung");
			sview.set_title(title);
			sview.show(u_new);


ogProjection.project_global(space,u_new, coeff_vec,  HERMES_L2_NORM); 
			Solution<double>::vector_to_solution(coeff_vec, space, u_taylor);
pview.show(u_taylor);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

