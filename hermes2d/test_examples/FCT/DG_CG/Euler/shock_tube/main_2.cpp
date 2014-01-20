#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "interface.h"
#include "euler_flux.h"





using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;


const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const double time_step = 5e-4;
const double T_FINAL = 0.231;                       // Time interval length.                   
// Equation parameters.  
 
     
// Kappa.
const double KAPPA = 1.4; 
// Inlet x-velocity (dimensionless).
const double V1_EXT =0.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

bool serendipity = true;

SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));

/*
SpaceSharedPtr<double> space_rho(new L2_NEW_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_NEW_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_NEW_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_NEW_Space<double>(mesh, P_INIT, serendipity));
*/
/*BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
 bview.show(space_rho);
  View::wait();*/

/*	SpaceSharedPtr<double> space_rho(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new SpaceBB<double>(mesh, P_INIT));*/


	int dof_rho = space_rho->get_num_dofs();
	int dof_v_x = space_rho_v_x->get_num_dofs();
	int dof_v_y = space_rho_v_y->get_num_dofs();
	int dof_e = space_e->get_num_dofs();

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));	


  MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_e(new Solution<double>);


  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> boundary_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);

double current_time =0;
 //--------- Visualization of pressure & velocity
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
 ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
 ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
//ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  ScalarView mach_view("mach", new WinGeom(0, 700, 600, 300));
			/*s2.set_min_max_range(1.,4.);
			s3.set_min_max_range(-1.,1.);
			mach_view.set_min_max_range(2.5,3.3);
			pressure_view.set_min_max_range(0.,3.);
			s1.set_min_max_range(1., 2.);
*/
	//OrderView m1view("mesh", new WinGeom(1000, 0, 500, 400));m1view.show(spaces[0]);

/*
MeshFunctionSharedPtr<double> mach_init(new  MachNumberFilter(init_slns, KAPPA));
		s1.show(init_rho);
		mach_view.show(mach_init);
*/
//View::wait(HERMES_WAIT_KEYPRESS);

//------------
	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
 
	//NumericalFlux* num_flux = new HLLNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new ApproxRoeNumericalFlux(KAPPA, euler_fluxes); 
NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new StegerWarmingNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new VijayasundaramNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new OsherSolomonNumericalFlux(KAPPA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);
	EulerK wf_convection_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA,prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

	EulerS wf_bdry_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerS wf_bdry(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);



  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp(&wf_convection, spaces); 
  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 

  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double>*  matrix = new CSCMatrix<double>;  
	CSCMatrix<double>*  mat_rhs = new CSCMatrix<double>; 
	CSCMatrix<double>* mass_matrix = new CSCMatrix<double>;  
	CSCMatrix<double>* dS_matrix = new CSCMatrix<double>;  
CSCMatrix<double>* dg_matrix = new CSCMatrix<double>; 
CSCMatrix<double>* K_matrix = new CSCMatrix<double>; 
CSCMatrix<double>* K_matrix_new = new CSCMatrix<double>; 

CSCMatrix<double>*   K_n= new CSCMatrix<double>;
CSCMatrix<double>*   K_n_1= new CSCMatrix<double>; ;

    OGProjection<double> ogProjection;
		SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_bdry = new SimpleVector<double> (ndof);	

		double* coeff_vec = new double[ndof];	
		double* coeff_vec_2 = new double[ndof];
double* mc_u_n  = new double[ndof];
double* mc_u_n_1  = new double[ndof];

		SimpleVector<double> * s_n = new SimpleVector<double> (ndof);
		SimpleVector<double> * s_n_1 = new SimpleVector<double> (ndof);

//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	


MeshFunctionSharedPtr<double> mach_2(new  MachNumberFilter(prev_slns, KAPPA));
		//s1.show(prev_rho);
		//mach_view.show(mach_2);
//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	
	int ts = 1;
	char title[100];

double t_k = time_step/2.;
double t_k_1= time_step;	
double tau_1 = (t_k+t_k_1)/(t_k_1);
double tau_2 = t_k/t_k_1;

		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(mass_matrix);

//Timestep loop
do
{	 
  	//Hermes::Mixins::Loggable::Static::info("dofs=%i,%i,%i,%i= %i", dof_rho, dof_v_x,dof_v_y, dof_e, ndof);
  	Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof);
 	  
 	 if(ts!=1){
		for(int i = 0; i<ndof; i++) mc_u_n_1[i] = -t_k/(t_k*(t_k+t_k_1))*mc_u_n[i];
		mass_matrix->multiply_with_vector(coeff_vec, mc_u_n); 
		K_matrix->free();
		K_matrix->create(K_matrix_new->get_size(),K_matrix_new->get_nnz(), K_matrix_new->get_Ap(), K_matrix_new->get_Ai(),K_matrix_new->get_Ax());
		s_n_1->zero();
		s_n_1->add_vector(vec_bdry);
		//s_n_1->add_vector(vec_dg);
		dp.assemble(dg_matrix);
		dp_bdry.assemble(vec_bdry);
	  	dp_DG.assemble(K_matrix_new,vec_dg);
		K_matrix_new->add_sparse_matrix(dg_matrix);
		s_n->zero();
		s_n->add_vector(vec_bdry);
		//s_n->add_vector(vec_dg);
		s_n->multiply_with_Scalar(tau_1);
		s_n_1->multiply_with_Scalar(-tau_2);
	}else{
			dp_init.assemble(dg_matrix);
			dp_bdry_init.assemble(vec_bdry);
		 	dp_DG_init.assemble(K_matrix_new,vec_dg);	
			K_matrix_new->add_sparse_matrix(dg_matrix);
			K_matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());
			

			mass_matrix->multiply_with_vector(coeff_vec, mc_u_n);
			for(int i = 0; i<ndof; i++) mc_u_n_1[i] = -t_k/(t_k*(t_k+t_k_1))*mc_u_n[i]; 
			s_n->zero(); s_n_1->zero();
			s_n->add_vector(vec_bdry);
			//s_n->add_vector(vec_dg);
			s_n_1->add_vector(s_n);
			s_n->multiply_with_Scalar(tau_1);
			s_n_1->multiply_with_Scalar(-tau_2);
		}


		K_n_1->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());
		K_n_1->multiply_with_Scalar(tau_2); 
		K_n->create(K_matrix_new->get_size(),K_matrix_new->get_nnz(), K_matrix_new->get_Ap(), K_matrix_new->get_Ai(),K_matrix_new->get_Ax());
		K_n->multiply_with_Scalar(-tau_1); 


		mat_rhs->create(mass_matrix->get_size(),mass_matrix->get_nnz(), mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
		mat_rhs->multiply_with_Scalar((2.*t_k+t_k_1)/((t_k_1+t_k)*t_k)); 
		matrix->create(K_n->get_size(),K_n->get_nnz(), K_n->get_Ap(),K_n->get_Ai(),K_n->get_Ax());
		matrix->add_sparse_matrix(K_n_1);		
		matrix->add_sparse_matrix(mat_rhs);	


	//-------------rhs: M/tau+ (1-theta)(L) u^n------------	
		
			vec_rhs->zero(); 
			vec_rhs->add_vector(mc_u_n); 
			vec_rhs->multiply_with_Scalar((t_k+t_k_1)/(t_k*t_k_1));
			vec_rhs->add_vector(mc_u_n_1); 
			vec_rhs->add_vector(s_n); 
			vec_rhs->add_vector(s_n_1); 
		


	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}			

		for(int i=0; i<ndof;i++)
					coeff_vec[i] = solver->get_sln_vector()[i];		

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	

			// Visualize the solution.
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(prev_slns));
			sprintf(title, "vx: ts=%i",ts);	
			s2.set_title(title);
			sprintf(title, "vy: ts=%i",ts);	
			s3.set_title(title);
			vel_x->reinit();
			vel_y->reinit();
			s2.show(vel_x);
			s3.show(vel_y);
			MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
			sprintf(title, "Pressure: ts=%i",ts);
			pressure_view.set_title(title);
			pressure->reinit();
			pressure_view.show(pressure);

			MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));
			sprintf(title, "Mach: ts=%i",ts);
			mach_view.set_title(title);
			mach->reinit();
			mach_view.show(mach);

			sprintf(title, "Density: ts=%i",ts);
			s1.set_title(title);
			s1.show(prev_rho);



	//View::wait(HERMES_WAIT_KEYPRESS);

	  // Update global time.
  current_time += time_step;


  // Increase time step counter
  ts++;

		delete solver;
		matrix->free();
		mat_rhs->free();
		K_n->free();
		K_n_1->free();

}while (current_time < T_FINAL);




		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));

				pressure->reinit();
				mach->reinit();


        Linearizer lin_p;
			lin_p.save_solution_vtk(pressure, "p_end.vtk", "pressure", true);
        Linearizer lin_m;
			lin_m.save_solution_vtk(mach, "m_end.vtk", "mach", true);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_slns[0], "rho_end.vtk", "density", true);



Orderizer ord_space;
ord_space.save_orders_vtk(spaces[0], "space.vtk");



		//Cleanup
			delete mass_matrix;
			delete matrix; 
			delete mat_rhs;
			delete dS_matrix;  
			delete dg_matrix;
			delete K_matrix;
			delete K_matrix_new;
			delete K_n; 
			delete K_n_1;

			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;
			delete vec_bdry;



			delete num_flux;






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

