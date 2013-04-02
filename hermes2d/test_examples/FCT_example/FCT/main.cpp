#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include <list>

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =7;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const int P_MAX = 2; 
const double h_max = 0.1;                       
const double time_step = 5e-4;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length. 
 

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";



int main(int argc, char* argv[])
{
   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) 
  		basemesh.refine_all_elements();
  mesh.copy(&basemesh);


  // Initialize boundary conditions.
  //DefaultEssentialBCConst<double>  bc_essential(BDY_IN, 0.0);
  //EssentialBCs<double>  bcs(&bc_essential);
  
  // Create an space with default shapeset.  
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


  
    // Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;
	char title[100];
	

	
	int ref_ndof = space.get_num_dofs();	
	info(" ndof = %d ", ref_ndof); 
  
        // Initialize matrix solver.
     UMFPackMatrix<double> * matrix = new UMFPackMatrix<double> ;  
		UMFPackVector<double> * rhs = new UMFPackVector<double> (ref_ndof); 
    UMFPackLinearSolver<double>* solver = new UMFPackLinearSolver<double>( matrix, rhs);  


    // Initialize the FE problem.    
 	 CustomWeakForm wf(time_step, theta, &u_prev_time);
    DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, &space);  
    
     dp->assemble(matrix, rhs);      
     if(solver->solve()) Solution<double>::vector_to_solution(solver->get_sln_vector(), &space, &u_new);
    else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		//sview.show(&u_new);    
		//View::wait(HERMES_WAIT_KEYPRESS);

    CustomWeakForm wf_2(time_step, theta, &u_new);  
        double current_time = time_step;
        // Initialize the FE problem.
    DiscreteProblem<double>* dp_2 = new DiscreteProblem<double>(&wf_2, &space);
    int ts =1;
  do
  {     
  info(" Time step %d, time %3.5f", ts, current_time); 	
     dp_2->assemble(matrix, rhs);      
     if(solver->solve()) Solution<double>::vector_to_solution(solver->get_sln_vector(), &space, &u_new);
		else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		
	/*	sprintf(title, "Loesung: Time %3.2f,timestep %i", current_time,ts);
			 sview.set_title(title);
		sview.show(&u_new);*/
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
	 	
	 	
	 	
	 	
	 	
	  }
  while (current_time<T_FINAL);
  
 
  	//	sview.show(&u_new);
CustomInitialCondition exact_solution(space.get_mesh());
	Adapt<double>* error_estimation = new Adapt<double>(&space, HERMES_L2_NORM);	
double err_est = error_estimation->calc_err_est(&exact_solution,&u_new,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
double err_est_2 = error_estimation->calc_err_est(&u_new,&exact_solution,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);

info("err_est = %f, err_est_2 =%f,  ndof = %d", err_est,err_est_2, ref_ndof);

FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "err_est = %f, err_est_2 =%f,  ndof = %d", err_est,err_est_2, ref_ndof);
fclose (pFile);  
  
  
  
  
  
  

	delete dp;
	delete dp_2;
	delete matrix;
	delete rhs;
	delete solver;




  // Wait for the view to be closed.
  View::wait();
  return 0;
}

