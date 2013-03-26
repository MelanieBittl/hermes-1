#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =7;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.

                    
const double time_step = 5e-4;                           // Time step.
const double T_FINAL = 2.5;  

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

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
  mloader.load("domain_2.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  
  // Create an space with default shapeset.  

//L2_SEMI_CG_Space<double> space(&mesh,P_INIT);
//L2Space<double> space(&mesh, P_INIT);
H1Space<double> space(&mesh, P_INIT);

  int ndof = space.get_num_dofs();
  
  /* BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
  bview.show(&space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh); 
 Solution<double> u_new;
 
  // Initialize views.
	ScalarView sview("Loesung", new WinGeom(500, 500, 500, 400));
	      sview.set_min_max_range(-0.01, 1.);
	ScalarView lview("anfangs-Loesung", new WinGeom(500, 0, 500, 400));
lview.show(&u_prev_time);
  OGProjection<double> ogProjection;
  
    // Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;
	char title[100];
    Hermes::HermesCommonApi.set_integral_param_value(Hermes::exceptionsPrintCallstack, 0);
    Hermes::Hermes2D::Hermes2DApi.set_integral_param_value(Hermes::Hermes2D::numThreads, 1);	
	
	int ref_ndof = space.get_num_dofs();	
	//Hermes::Mixins::Loggable::Static::info(" ndof = %d ", ref_ndof); 
	
	
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   //M_c/tau
	UMFPackMatrix<double> * conv_matrix = new UMFPackMatrix<double> ;   //K
	UMFPackMatrix<double>* dg_surface_matrix = new UMFPackMatrix<double> ; //inner and outer edge integrals
	UMFPackVector<double> * surf_rhs = new UMFPackVector<double> (ref_ndof); 
		UMFPackVector<double> * mass_rhs = new UMFPackVector<double> (ref_ndof); 
	CustomWeakFormMassmatrix  massmatrix(time_step, &u_prev_time);
	CustomWeakFormConvection  convection(&u_prev_time);
	///////////////////////////------------false, false (only CG), false, true (DG) --------------------------------------------------------------------------------
	CustomWeakForm wf_surf(time_step, theta, &u_prev_time, BDY_IN, &mesh,false, false);
		///////////////////////////--------------------------------------------------------------------------------------------
	DiscreteProblem<double> * dp_mass = new DiscreteProblem<double> (&massmatrix, &space);
	DiscreteProblem<double> * dp_convection = new DiscreteProblem<double> (&convection, &space);
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,&space);	
	dp_mass->assemble(mass_matrix, mass_rhs); 
	dp_convection->assemble(conv_matrix, NULL,true);
	dp_surf->assemble(dg_surface_matrix,surf_rhs);
	
		UMFPackMatrix<double> * high_rhs = new UMFPackMatrix<double> ; 
		UMFPackMatrix<double> * matrix = new UMFPackMatrix<double> ;
		     
			matrix->create(dg_surface_matrix->get_size(),dg_surface_matrix->get_nnz(), dg_surface_matrix->get_Ap(), dg_surface_matrix->get_Ai(),dg_surface_matrix->get_Ax());
			matrix->multiply_with_Scalar(-1.); //damit surface richtiges Vorzeichen!!
			matrix->add_matrix(conv_matrix); 
			matrix->multiply_with_Scalar(-theta);
			matrix->add_matrix(mass_matrix);  
			
			high_rhs->create(dg_surface_matrix->get_size(),dg_surface_matrix->get_nnz(), dg_surface_matrix->get_Ap(), dg_surface_matrix->get_Ai(),dg_surface_matrix->get_Ax());
			high_rhs->multiply_with_Scalar(-1.); //damit surface richtiges Vorzeichen!!
			high_rhs->add_matrix(conv_matrix); 
			high_rhs->multiply_with_Scalar((1.0-theta));
			high_rhs->add_matrix(mass_matrix); 
	
	double *coeff_vec = new double[ref_ndof];		
	double *coeff_vec_2 = new double[ref_ndof];
	double* vec_new;
  memset(coeff_vec_2, 0, ref_ndof*sizeof(double));

  
        // Initialize matrix solver.  
		UMFPackVector<double> * rhs = new UMFPackVector<double> (ref_ndof); 
		

    int ts = 1;
    double current_time =0.;

  do
  {     
  Hermes::Mixins::Loggable::Static::info(" Time step %d, time %3.5f", ts, current_time); 
  
     high_rhs->multiply_with_vector(coeff_vec_2, coeff_vec); 
		rhs->zero(); rhs->add_vector(surf_rhs);
		rhs->add_vector(coeff_vec);
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>( matrix, rhs); 
    
    
     if(solver->solve())
     {
     	vec_new = solver->get_sln_vector();
      Solution<double>::vector_to_solution(vec_new, &space, &u_new);
      }
    else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		
		for(int i = 0; i<ref_ndof; i++) 
			coeff_vec_2[i] = vec_new[i];
		
			sprintf(title, "Loesung: Time %3.2f,timestep %i", current_time,ts);
			 sview.set_title(title);
		sview.show(&u_new);
		//View::wait(HERMES_WAIT_KEYPRESS);	
	
	 	current_time += time_step;
	 	ts++;
	 	
	 	
	 	  // Visualization.
  /*  if((ts - 1) % 1000 == 0) 
    {
      // Output solution in VTK format.
        char filename[40];
        sprintf(filename, "solution-%i.vtk", ts );
        lin.save_solution_vtk(&u_new, filename, "solution", mode_3D);  
     
    }*/	
	 	
	 	delete solver;
	 	
	  }
  while (current_time<T_FINAL);
  

//lin.save_solution_vtk(&u_new, "sln.vtk", "solution", mode_3D);
  

CustomInitialCondition exact_solution(space.get_mesh());

/*
double bdry_err = calc_error_bdry(&u_new,&exact_solution, &space);
			sprintf(title, "Loesung: Time %3.2f,timestep %i", current_time,ts);
			 sview.set_title(title);
		sview.show(&u_new);
		//View::wait(HERMES_WAIT
*/

  ogProjection.project_global(&space, &exact_solution, coeff_vec,  HERMES_L2_NORM);  
 double abs_err_l2 = Global<double>::calc_abs_error(&exact_solution,&u_new, HERMES_L2_NORM);
 double err_l2 = calc_error_l2(&exact_solution, &u_new, &space);
 double abs_err_max = calc_error_max(&exact_solution, &u_new, &space);
 double vec_err_max = calc_error_max(coeff_vec, coeff_vec_2 ,ref_ndof);

Hermes::Mixins::Loggable::Static::info("l2=%f, l2_2= %f abs_max = %f, vec_max = %f, ndof = %d", abs_err_l2,err_l2, abs_err_max ,vec_err_max, ref_ndof);

FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%f,abs_max = %f, vec_max = %f, ndof = %d", abs_err_l2, abs_err_max, vec_err_max,  ref_ndof);
fclose (pFile);  
  
  
  
  
    delete [] coeff_vec;
  delete [] coeff_vec_2;


	delete matrix;
	delete rhs;
	
	delete mass_matrix;  
	delete conv_matrix;
	delete high_rhs;
	delete dg_surface_matrix;
		delete dp_convection;
		delete dp_mass; 
delete dp_surf;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

