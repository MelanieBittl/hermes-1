#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "interface.h"
#include "lumped_projection.h"
#include "euler_flux.h"
#include "hp_adapt.h"




using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;



const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const double time_step = 5e-4;
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


//FCT 
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
  mloader.load("domain.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));	


SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));*/

bool serendipity = true;

SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));


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


 /*   MeshFunctionSharedPtr<double> prev_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> prev_rho_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> prev_rho_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> prev_e(new CustomInitialCondition_e(mesh,KAPPA));	*/

  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> boundary_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	

		MeshFunctionSharedPtr<double> low_rho(new Solution<double>);
		MeshFunctionSharedPtr<double> low_rho_v_x(new Solution<double>);
		MeshFunctionSharedPtr<double> low_rho_v_y(new Solution<double>);
		MeshFunctionSharedPtr<double>	low_e(new Solution<double>);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > low_slns(low_rho,low_rho_v_x,low_rho_v_y,low_e);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));
      s1.set_min_max_range(0, 1.);
      s2.set_min_max_range(0., 1.);
      s3.set_min_max_range(0., 1.);
			pressure_view.set_min_max_range(0.,1.);

	OrderView m1view("mesh", new WinGeom(1000, 0, 500, 400));
	OrderView m2view("mesh", new WinGeom(1200, 0, 500, 400));
	OrderView m3view("mesh", new WinGeom(1400, 0, 500, 400));
	OrderView m4view("mesh", new WinGeom(1500, 0, 500, 400));

//------------
	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
 
	//NumericalFlux* num_flux = new HLLNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new ApproxRoeNumericalFlux(KAPPA, euler_fluxes); 
NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new StegerWarmingNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new VijayasundaramNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new OsherSolomonNumericalFlux(KAPPA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA,mesh, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA,mesh, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);

  EulerK   wf_K_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e);
  //EulerBoundary wf_boundary_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);


  EulerEquationsWeakForm_Mass wf_mass;
  EulerK   wf_K(KAPPA,  prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  //EulerBoundary wf_boundary(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);


  EulerK   wf_K_low(KAPPA,low_rho, low_rho_v_x, low_rho_v_y, low_e  );
 // EulerBoundary wf_boundary_low(KAPPA,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e,low_rho, low_rho_v_x, low_rho_v_y, low_e  );


EulerS wf_boundary_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);
EulerS wf_boundary(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
EulerS wf_boundary_low(KAPPA,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e,low_rho, low_rho_v_x, low_rho_v_y, low_e  );


  // Initialize the FE problem.
  DiscreteProblem<double> dp_boundary_init(&wf_boundary_init,spaces);
  DiscreteProblem<double> dp_K_init(&wf_K_init, spaces);
    
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_boundary(&wf_boundary, spaces);
  DiscreteProblem<double> dp_K(&wf_K, spaces);

  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);

  DiscreteProblem<double> dp_boundary_low(&wf_boundary_low, spaces);
  DiscreteProblem<double> dp_K_low(&wf_K_low, spaces);


  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> * matrix_L_low = new CSCMatrix<double>; 
	CSCMatrix<double> * low_matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * lowmat_rhs = new CSCMatrix<double>; 
	CSCMatrix<double> * matrix_dS = new CSCMatrix<double>; 
	CSCMatrix<double> * matrix_dS_low = new CSCMatrix<double>; 
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double>;  

	CSCMatrix<double> * ks_matrix = new CSCMatrix<double>;  

    OGProjection<double> ogProjection;
	Lumped_Projection lumpedProjection;

//Projection of the initial condition

	double* coeff_vec_init = new double[ndof];
  ogProjection.project_global(spaces,init_slns, coeff_vec_init, norms_l2 );
	Solution<double>::vector_to_solutions(coeff_vec_init, spaces, low_slns);	

		KrivodonovaDiscontinuityDetector dis_detect_init(spaces, low_slns);
		std::set<int> discont_elem =	dis_detect_init.get_discontinuous_element_ids();
		DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
		HPAdapt adapting(spaces, &error_calculator);
		bool new_adap = adapting.reduce_order(&discont_elem);


	//	m1view.show(spaces[0]);
		//m2view.show(spaces[1]);m3view.show(spaces[2]);m4view.show(spaces[3]);

	 delete [] coeff_vec_init; 



// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	

	double* u_L = NULL; 

int  ps =1; 
int ndof_alt = Space<double>::get_num_dofs(spaces);
		KrivodonovaDiscontinuityDetector* dis_detect = new KrivodonovaDiscontinuityDetector(spaces, low_slns);

			SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof_alt);
			SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof_alt);
			double* coeff_vec = new double[ndof_alt];	double* coeff_vec_2 = new double[ndof_alt];
			double* P_plus = new double[ndof_alt]; double* P_minus = new double[ndof_alt];
			double* Q_plus = new double[ndof_alt]; double* Q_minus = new double[ndof_alt];	
			double* R_plus = new double[ndof_alt]; double* R_minus = new double[ndof_alt];	

//Timestep loop
do
{	 

			ndof = Space<double>::get_num_dofs(spaces);
			bool* fct = get_fct_dofs(spaces);

			dof_rho = space_rho->get_num_dofs();
			dof_v_x = space_rho_v_x->get_num_dofs();
			dof_v_y = space_rho_v_y->get_num_dofs();
			dof_e = space_e->get_num_dofs();

if(ndof!=ndof_alt)
{
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;

			 vec_dg = new SimpleVector<double> (ndof);
			 vec_rhs = new SimpleVector<double> (ndof);
			 coeff_vec = new double[ndof];	coeff_vec_2 = new double[ndof];
			 P_plus = new double[ndof];  P_minus = new double[ndof];
			 Q_plus = new double[ndof];  Q_minus = new double[ndof];	
			 R_plus = new double[ndof];  R_minus = new double[ndof];	
			ndof_alt = ndof;

}


  	//Hermes::Mixins::Loggable::Static::info("dofs=%i,%i,%i,%i= %i", dof_rho, dof_v_x,dof_v_y, dof_e, ndof);
  	Hermes::Mixins::Loggable::Static::info("Time step %d,ps = %i,  time %3.5f, ndofs=%i", ts, ps,current_time, ndof);
		Space<double>::assign_dofs(spaces);	
   if((new_adap==true)||(ts==1))
{
			dp_mass.set_spaces(spaces);dp_boundary.set_spaces(spaces);dp_DG.set_spaces(spaces);dp_boundary_low.set_spaces(spaces);
		  dp_mass.assemble(mass_matrix);
}
		
		CSCMatrix<double> * lumped_matrix = massLumping(fct,mass_matrix);


 	  
 	  if(ts!=1)
 	  {
		//	dp_boundary.assemble(matrix_dS);
		 // dp_K.assemble(lowmat_rhs);
		dp_boundary.assemble(ks_matrix);

		  dp_DG.assemble(vec_dg);	

			lumpedProjection.project_lumped(spaces,prev_slns, coeff_vec, matrix_solver, fct);
			ogProjection.project_global(spaces, prev_slns, coeff_vec_2, norms_l2);

		}else{
		//	dp_boundary_init.assemble(matrix_dS);
		//  dp_K_init.assemble(lowmat_rhs);	
		dp_boundary_init.assemble(ks_matrix);

		  dp_DG_init.assemble(vec_dg);	
				lumpedProjection.project_lumped(spaces,init_slns, coeff_vec, matrix_solver, fct);
				ogProjection.project_global(spaces,init_slns, coeff_vec_2, norms_l2);

		}


		lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e,fct);
	/*Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);
  MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
 MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(prev_slns));	
			pressure->reinit();
			vel_x->reinit();
			vel_y->reinit();
			s1.show(prev_rho);     
			s2.show(vel_x);
			s3.show(vel_y);
  		pressure_view.show(pressure);
View::wait(HERMES_WAIT_KEYPRESS);*/

	//------------------------artificial DIFFUSION D---------------------------------------	
		Space<double>::assign_dofs(spaces);	
			CSCMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,spaces,ks_matrix, fct);
			ks_matrix->add_sparse_matrix(diffusion); //L(U)=K+D
			//lowmat_rhs->add_sparse_matrix(matrix_dS); //L(U)+dS(U) 
      lowmat_rhs->create_merged_pattern(ks_matrix,lumped_matrix);
			lowmat_rhs->add_sparse_matrix(ks_matrix);


			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_sparse_matrix(lumped_matrix); 				//M_L/t - theta L(U)
			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_sparse_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)


	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 
			vec_dg->multiply_with_Scalar(time_step);
			vec_rhs->add_vector(vec_dg); 

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
			try{
			 lowOrd->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}		
			u_L = lowOrd->get_sln_vector();  

      Solution<double>::vector_to_solutions(u_L, spaces,low_slns);
ps++;
if(ps==1){
		std::set<int> discont_elem =	dis_detect->get_discontinuous_element_ids();
		new_adap = adapting.reduce_order(&discont_elem);

	/*	if(new_adap==true) 
		{
				m1view.show(spaces[0]);	m2view.show(spaces[1]);m3view.show(spaces[2]);m4view.show(spaces[3]);		

			View::wait(HERMES_WAIT_KEYPRESS);
		}*/
}
if((ps==2)||(new_adap ==false))
{
			//dp_boundary_low.assemble(matrix_dS_low);	
    	//dp_K_low.assemble(matrix_L_low);

				dp_boundary_low.assemble(matrix_L_low);	

			CSCMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,spaces,matrix_L_low,fct);
			matrix_L_low->add_sparse_matrix(diffusion_low); //L(U)
			//matrix_L_low->add_sparse_matrix(matrix_dS_low); //L(U)+dS(U) 

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		//Hermes::Mixins::Loggable::Static::info("antidiffusive fluxes ");
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e, fct);

			lumped_matrix->multiply_with_vector(u_L, coeff_vec); 

		for(int i=0; i<ndof;i++)
					coeff_vec[i] += coeff_vec_2[i]*time_step;	

			vec_rhs->zero();
			vec_rhs->add_vector(coeff_vec); 
			UMFPackLinearMatrixSolver<double> * newsol = new UMFPackLinearMatrixSolver<double> (lumped_matrix,vec_rhs);	
			try{
			 newsol->solve();
			}catch(Hermes::Exceptions::Exception e)	{
				e.print_msg();
			}	

		for(int i=0; i<ndof;i++)
					coeff_vec[i] = newsol->get_sln_vector()[i];


			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	

			// Visualize the solution.

		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(prev_slns));	
			sprintf(title, " ts=%i",ts);
			pressure_view.set_title(title);
			s1.set_title(title);
			s2.set_title(title);
			s3.set_title(title);
			pressure->reinit();
			vel_x->reinit();
			vel_y->reinit();
			s1.show(prev_rho);
			s2.show(vel_x);
			s3.show(vel_y);
			pressure_view.show(pressure);
		m1view.show(spaces[0]);
	//View::wait(HERMES_WAIT_KEYPRESS);


		delete diffusion_low;
		delete newsol;

		ps =1.;
	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;


 }else ps++; 	


		delete lowOrd;
		delete diffusion;	
		delete lumped_matrix;
		low_matrix->free();
		lowmat_rhs->free();
		delete [] fct;

}
while (current_time < T_FINAL);




		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(prev_slns));	

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



Orderizer ord_space;
ord_space.save_orders_vtk(spaces[0], "space.vtk");



		//Cleanup
			delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete low_matrix;
			delete lowmat_rhs;
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;
delete dis_detect;





  // Wait for the view to be closed.
  View::wait();
  return 0;
}

