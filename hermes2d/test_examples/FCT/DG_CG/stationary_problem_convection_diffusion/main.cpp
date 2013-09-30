#define HERMES_REPORT_ALL
#include "definitions.h"
#include "hermes2d.h"
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.

///triangles : ref 6 => h = 1/64

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY = "bdry";



bool all = true;
bool DG = true;
bool SD = false;


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


EssentialBCNonConst bc_essential(BDY, u_prev_time);
EssentialBCs<double> bcs(&bc_essential);


  // Create an space with default shapeset.
  //SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT, serendipity));	
 SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh,P_INIT));

 /* BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
  bview.show(space);
View::wait(HERMES_WAIT_KEYPRESS);*/





  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	//lview.show(u_prev_time);
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
/*	OGProjection<double> ogProjection;	
				ogProjection.project_global(space_h1, u_prev_time, coeff_vec, HERMES_L2_NORM);		
  MeshFunctionSharedPtr<double>  u_test(new Solution<double>);
Solution<double>::vector_to_solution(coeff_vec, space_h1, u_test);*/


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
			for(int i=0; i<ndof; i++) coeff_vec_2[i] = vec_new[i];

		sview.show(u_new);

/*
AsmList<double>* al  = new AsmList<double>;		Element* e;
 MeshFunctionSharedPtr<double>  u_test(new Solution<double>);
	for_all_active_elements(e, space->get_mesh()) 
	{
				space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		 			 int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_TRIANGLE);
								Node* vn = e->vn[iv];
								if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
									for(unsigned int j = 0; j < al->get_cnt(); j ++){	
											int dof = al->get_dof()[j];		 
											if((al->get_idx()[j]==index)&&(dof!=-1.0)){ 
													coeff_vec_3[dof]= coeff_vec_2[dof];
											}
									}
							 	}
							}	


		}

			Solution<double>::vector_to_solution(coeff_vec_3, space, u_test);
	lview.show(u_test);*/

calc_error_total(u_new, u_prev_time,space);






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

