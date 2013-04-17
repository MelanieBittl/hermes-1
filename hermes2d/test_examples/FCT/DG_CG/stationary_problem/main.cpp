#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.

                    
const double time_step = 1.;                           // Time step.
const double T_FINAL = 3.0;  

const double theta = 1.;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

//FCT & p-Adaptivity
#include "error_estimations.cpp"

int main(int argc, char* argv[])
{
   // Load the mesh->
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  
  // Create an space with default shapeset.  

//CustomDirichletCondition bc_essential(Hermes::vector<std::string>("inlet1","inlet2"));
 // EssentialBCs<double>  bcs(&bc_essential);

L2_SEMI_CG_Space<double> space(&mesh,P_INIT);
//L2Space<double> space(&mesh, P_INIT);
//H1Space<double> space(&mesh,&bcs, P_INIT);
//H1Space<double> space(&mesh,P_INIT);

  int ndof = space.get_num_dofs();
  
 /*  BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
  bview.show(&space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh); 
 Solution<double> u_new, proj_sln;
 
  // Initialize views.
	ScalarView sview("Loesung", new WinGeom(500, 500, 500, 400));
ScalarView fview("filter", new WinGeom(500, 500, 500, 400));
	      //sview.set_min_max_range(-0.01, 1.);
	ScalarView lview("anfangs-Loesung", new WinGeom(500, 0, 500, 400));
//lview.show(&u_prev_time);
  OGProjection<double> ogProjection;
  
    // Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;
	char title[100];
    Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
    Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);	
	
	int ref_ndof = space.get_num_dofs();	
	//Hermes::Mixins::Loggable::Static::info(" ndof = %d ", ref_ndof); 


	double *coeff_vec = new double[ref_ndof];		
	double *coeff_vec_2 = new double[ref_ndof];
	double *coeff_vec_3 = new double[ref_ndof];
	double* vec_new;
  memset(coeff_vec_2, 0, ref_ndof*sizeof(double));

CustomWeakFormConvection  convection;
	///////////////////////////------------false, false (only CG), false, true (DG) --------------------------------------------------------------------------------
	CustomWeakForm wf_surf(time_step, theta, &u_prev_time, BDY_IN, &mesh,true,true);
		///////////////////////////--------------------------------------------------------------------------------------------


	UMFPackMatrix<double>* dg_surface_matrix = new UMFPackMatrix<double> ; //inner and outer edge integrals
	UMFPackVector<double> * surf_rhs = new UMFPackVector<double> (ref_ndof); 
	DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, &space);
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,&space);	
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	dp_convection->assemble(conv_matrix, NULL,true);

	dp_surf->assemble(dg_surface_matrix,surf_rhs);
//dp_surf->assemble(dg_surface_matrix);
	
		UMFPackMatrix<double> * high_rhs = new UMFPackMatrix<double> ; 
		UMFPackMatrix<double> * matrix = new UMFPackMatrix<double> ;
		     
			matrix->create(dg_surface_matrix->get_size(),dg_surface_matrix->get_nnz(), dg_surface_matrix->get_Ap(), dg_surface_matrix->get_Ai(),dg_surface_matrix->get_Ax());
			//matrix->add_matrix(conv_matrix); 
  
        // Initialize matrix solver.  
		//UMFPackVector<double> * rhs = new UMFPackVector<double> (ref_ndof);   
 
		//rhs->zero(); rhs->add_vector(surf_rhs);
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>( matrix, surf_rhs); 
    
    
     if(solver->solve())
     {
     	vec_new = solver->get_sln_vector();
      Solution<double>::vector_to_solution(vec_new, &space, &u_new);
      }
    else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");

		for(int i = 0; i<ref_ndof; i++) 
			coeff_vec_2[i] = vec_new[i];

		sview.show(&u_new);
		//View::wait(HERMES_WAIT_KEYPRESS);	

/*	 	matrix->multiply_with_vector(coeff_vec_2, coeff_vec); 
	for(int i =0; i<ref_ndof; i++) coeff_vec_3[i] = coeff_vec[i] - surf_rhs->get(i);
double total = 0;
	for(int i =0;i<ref_ndof; i++) total += abs(coeff_vec_3[i]);	 	

Hermes::Mixins::Loggable::Static::info("total = %.5e      ", total);
	 	  // Visualization.
  /*  if((ts - 1) % 1000 == 0) 
    {
      // Output solution in VTK format.
        char filename[40];
        sprintf(filename, "solution-%i.vtk", ts );
        lin.save_solution_vtk(&u_new, filename, "solution", mode_3D);  
     
    }*/	
	 	
	 	//delete solver;
	 	
 // Initialize linear solver.
/*  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf_surf, &space);

  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, &space, &u_new);

		for(int i = 0; i<ref_ndof; i++) 
			coeff_vec_2[i] = sln_vector[i];

		sview.show(&u_new);
  }catch(std::exception& e)
  {
    std::cout << e.what();
  }
 */ 


 

CustomInitialCondition exact_solution(space.get_mesh());

  ogProjection.project_global(&space, &exact_solution, coeff_vec,  HERMES_L2_NORM);  
Solution<double>::vector_to_solution(coeff_vec, &space, &proj_sln);

 double abs_err_l2 = Global<double>::calc_abs_error(&exact_solution,&u_new, HERMES_L2_NORM);
 double err_l2 = calc_error_l2(&proj_sln, &u_new, &space);
 double err_l1 = calc_error_l1(&exact_solution, &u_new, &space);
 double err_l1_proj = calc_error_l1(&proj_sln, &u_new, &space);
 double abs_err_max = calc_error_max(&exact_solution, &u_new, &space);
 double err_max_coeff = calc_error_max(coeff_vec, coeff_vec_2, ref_ndof);
Hermes::Mixins::Loggable::Static::info("l2=%.5e, l2_new = %.5e, l1=%.5e, abs_max = %f,  ndof = %d", abs_err_l2,err_l2,err_l1, abs_err_max , ref_ndof);

AbsDifffilter filter(Hermes::vector<MeshFunction<double>*>(&exact_solution, &u_new));
//fview.show(&filter);
lin.save_solution_vtk(&u_new, "sln.vtk", "solution", mode_3D);
lin.save_solution_vtk(&filter, "error.vtk" , "error", mode_3D);  

FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%.5e, l2(proj) = %.5e, l1=%.5e,l1(proj)=%.5e, abs_max = %.5e,abs_coeff = %.5e  ndof = %d", abs_err_l2,err_l2,err_l1,err_l1_proj, abs_err_max , err_max_coeff, ref_ndof);
fclose (pFile);  

 
  
  
    delete [] coeff_vec;
  delete [] coeff_vec_2;


/*	delete matrix;
	delete rhs;

	delete high_rhs;
	delete dg_surface_matrix;
	delete dp_surf;*/


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

