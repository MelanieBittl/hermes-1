#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "prev_solution.h"
#include"util.h"
#include "shapeset_taylor.h" 
 #include "solution_slopelimiter.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 


const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const int P_MAX = 2; 
const double h_max = 0.1;                       
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.
 //const double T_FINAL = 1e-3;    

 
const double EPS = 1e-12;
const double EPS_smooth = 1e-14;

const double THRESHOLD = 0.3;

const int NDOF_STOP = 20000;   

const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

//FCT & p-Adaptivity
#include "fct.cpp"
#include "p_adapt.cpp"
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "p1_list.cpp"
#include "reg_estimator.cpp"
#include "z_z.cpp"



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

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));	 


  // Initialize views.
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));




  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
Orderizer ord;




// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int as = 1;  //p-adapt-schritte


	 double diag;			double f;
	Element* e =NULL;
	for_all_active_elements(e, mesh){diag = e->get_diameter(); break;}
	const double h_start = diag;
	const double h_min = diag/8; 


	AsmList<double>*  al = new AsmList<double>;	

  OGProjection<double> ogProjection;

	//SpaceSharedPtr<double> ref_space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
		
SpaceSharedPtr<double> ref_space_2(new L2Space<double>(mesh,P_INIT));	

int ref_ndof = ref_space_2->get_num_dofs(); 	
   double* coeff_vec = new double[ref_ndof];
		ogProjection.project_global(ref_space_2,u_prev_time, coeff_vec,  HERMES_L2_NORM); 
			Solution<double>::vector_to_solution(coeff_vec, ref_space_2, u_new);
	

			//sview.show(u_prev_time);

SlopeLimiterSolution sln(u_new, ref_space_2);
sln.limit_solution_according_to_detector();


		sprintf(title, "lim. Loesung");
			pview.set_title(title);
			pview.show(&sln);






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

