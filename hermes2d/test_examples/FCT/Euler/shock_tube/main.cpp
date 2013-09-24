#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "lumped_projection.h"
#include <list>

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const double time_step = 125e-6;
const double T_FINAL = 0.231;                       // Time interval length. 

const double theta = 0.5;

// Equation parameters.  
 
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;  

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"    

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain2.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

	/*	SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));	*/


		SpaceSharedPtr<double> space_rho(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new SpaceBB<double>(mesh, P_INIT));

		int dof_rho = space_rho->get_num_dofs();
		int dof_v_x = space_rho_v_x->get_num_dofs();
		int dof_v_y = space_rho_v_y->get_num_dofs();
		int dof_e = space_e->get_num_dofs();

    Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
		Space<double>::assign_dofs(spaces);
  	int ndof = Space<double>::get_num_dofs(spaces);
  	printf("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
		MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh));	
		MeshFunctionSharedPtr<double> init_rho_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
		MeshFunctionSharedPtr<double> init_rho_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
		MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));	

		MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
		MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
		MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
		MeshFunctionSharedPtr<double> prev_e(new Solution<double>);

		MeshFunctionSharedPtr<double> low_rho(new Solution<double>);
		MeshFunctionSharedPtr<double> low_rho_v_x(new Solution<double>);
		MeshFunctionSharedPtr<double> low_rho_v_y(new Solution<double>);
		MeshFunctionSharedPtr<double>	low_e(new Solution<double>);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > low_slns(low_rho,low_rho_v_x,low_rho_v_y,low_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);

 //--------- Visualization of pressure & velocity
/*  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));
      s1.set_min_max_range(0, 1.);
      s2.set_min_max_range(0., 1.);
      s3.set_min_max_range(0, 0.1);
			pressure_view.set_min_max_range(0.,1.);*/

//--------------Weakforms------------
  EulerEquationsWeakForm_K  wf_K_init(KAPPA, init_slns);
  EulerBoundary wf_boundary_init(KAPPA, Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e,init_rho, init_rho_v_x, init_rho_v_y,  init_e) );
  EulerEquationsWeakForm_Mass wf_mass;

  EulerEquationsWeakForm_K*  wf_K= new EulerEquationsWeakForm_K(KAPPA, prev_slns);
  EulerBoundary* wf_boundary= new EulerBoundary(KAPPA,
Hermes::vector<MeshFunctionSharedPtr<double> > (prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, init_rho, init_rho_v_x, init_rho_v_y,  init_e));
  EulerEquationsWeakForm_K*  wf_K_low = new EulerEquationsWeakForm_K(KAPPA, low_slns);
  EulerBoundary* wf_boundary_low = new EulerBoundary(KAPPA,Hermes::vector<MeshFunctionSharedPtr<double> > (low_rho, low_rho_v_x, low_rho_v_y, low_e, init_rho, init_rho_v_x, init_rho_v_y,  init_e ));


  //------------------- Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_boundary_init(&wf_boundary_init,spaces);
  DiscreteProblem<double> dp_K_init(&wf_K_init, spaces);
  DiscreteProblem<double>* dp_boundary = new DiscreteProblem<double>(wf_boundary, spaces);
  DiscreteProblem<double>* dp_K= new DiscreteProblem<double>(wf_K, spaces);
  DiscreteProblem<double>* dp_boundary_low= new DiscreteProblem<double>(wf_boundary_low, spaces);
  DiscreteProblem<double>* dp_K_low= new DiscreteProblem<double>(wf_K_low, spaces);


  // Set up the matrix according to the solver selection.
 
	CSCMatrix<double> * matrix_L_low = new CSCMatrix<double>; 
	CSCMatrix<double> * low_matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * lowmat_rhs = new CSCMatrix<double>; 
	CSCMatrix<double> * matrix_dS = new CSCMatrix<double>; 
	CSCMatrix<double> * matrix_dS_low = new CSCMatrix<double>; 
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double>;  
	SimpleVector<double> * vec_rhs = new SimpleVector<double>(ndof);
	double* u_L = NULL; 

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	



///-----------Assembling mass matrix
    dp_mass.assemble(mass_matrix);
	//----------------------MassLumping M_L--------------------------------------------------------------------
		CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);
//---------------------Projection of the initial condition
    OGProjection<double> ogProjection;
		Lumped_Projection lumpedProjection;

	lumpedProjection.project_lumped(spaces, init_slns, coeff_vec, matrix_solver);
  ogProjection.project_global(spaces,init_slns, coeff_vec_2, norms_l2 );
	lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);
	Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);



//------------Filter---------------------------
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double>  vel_x(new VelocityFilter(prev_slns, 1));
  MeshFunctionSharedPtr<double>  vel_y(new VelocityFilter(prev_slns, 2));

			 // Visualize the solution.
			//s1.show(prev_rho);
//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];


Space<double>::assign_dofs(spaces);
//Timestep loop
do
{	 
  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time); 	  
 	  if(ts!=1) {
			dp_boundary->assemble(matrix_dS);
		  dp_K->assemble(lowmat_rhs);
		}else{
			dp_boundary_init.assemble(matrix_dS);
		  dp_K_init.assemble(lowmat_rhs);				
		}
	//------------------------artificial DIFFUSION D---------------------------------------	
			CSCMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);

			lowmat_rhs->add_sparse_matrix(diffusion); //L(U)=K+D
			lowmat_rhs->add_sparse_matrix(matrix_dS); //L(U)+dS(U) 

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_sparse_matrix(lumped_matrix); 				//M_L/t - theta L(U)

			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_sparse_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)

	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
			try{
			 lowOrd->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}		
			u_L = lowOrd->get_sln_vector();  
      Solution<double>::vector_to_solutions(u_L, spaces,low_slns);

			dp_boundary_low->assemble(matrix_dS_low);	
    	dp_K_low->assemble(matrix_L_low);
			CSCMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_sparse_matrix(diffusion_low); //L(U)
			matrix_L_low->add_sparse_matrix(matrix_dS_low); //L(U)+dS(U) 

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		//Hermes::Mixins::Loggable::Static::info("antidiffusive fluxes ");
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

		for(int i=0; i<ndof;i++)
					coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);					

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	
	

			 // Visualize the solution.
		/*		
				pressure->reinit();
				vel_x->reinit();
				vel_y->reinit();
				sprintf(title, "pressure: ts=%i",ts);
				pressure_view.set_title(title);
				pressure_view.show(pressure);
				sprintf(title, "density: ts=%i",ts);
				s1.set_title(title);
				s1.show(prev_slns[0]);
				sprintf(title, "velocity (x): ts=%i",ts);
				s2.set_title(title);
				s2.show(vel_x);
				sprintf(title, "velocity (y): ts=%i",ts);
				s3.set_title(title);
				s3.show(vel_y);

  		

		//	View::wait(HERMES_WAIT_KEYPRESS);



 /* // Visualization.
    if((ts - 1) % EVERY_NTH_STEP == 0) 
    {

        pressure->reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-%i.vtk", ts - 1);
        lin_pressure.save_solution_vtk(pressure, filename, "Pressure", true);

      
    }*/




		delete lowOrd;
		delete diffusion;
		delete diffusion_low;
		low_matrix->free();



	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);

				pressure->reinit();
				vel_x->reinit();
				vel_y->reinit();
        Linearizer lin_p;
			lin_p.save_solution_vtk(pressure, "p_end.vtk", "pressure", true);
        Linearizer lin_v_x;
			lin_v_x.save_solution_vtk(vel_x, "vx_end.vtk", "velocity_x", true);
        Linearizer lin_v_y;
			lin_v_y.save_solution_vtk(vel_y, "vy_end.vtk", "velocity_y",true);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_slns[0], "rho_end.vtk", "density", true);


delete dp_boundary; delete dp_K; delete wf_K; delete wf_boundary;

delete dp_boundary_low;
delete wf_boundary_low;
delete dp_K_low;
delete wf_K_low;

		//Cleanup
		delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete lumped_matrix;
			delete low_matrix;
			delete lowmat_rhs;
			delete vec_rhs;
			  // Clean up.
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete[] coeff_vec_2;
			delete [] coeff_vec; 



  // Wait for the view to be closed.
  View::wait();
  return 0;
}

