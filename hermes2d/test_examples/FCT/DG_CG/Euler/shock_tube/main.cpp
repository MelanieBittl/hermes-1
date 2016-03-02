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



const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 2;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 0.231;                       // Time interval length. 

const double theta = 1.;
const bool DG = true;
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
  //mloader.load("domain2.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));
/*		
SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));

/*SpaceSharedPtr<double> spacel2_rho(new L2_NEW_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> spacel2_rho_v_x(new L2_NEW_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> spacel2_rho_v_y(new L2_NEW_Space<double>(mesh, P_INIT));	
SpaceSharedPtr<double> spacel2_e(new L2_NEW_Space<double>(mesh, P_INIT));
Hermes::vector<SpaceSharedPtr<double> > spacesl2(spacel2_rho, spacel2_rho_v_x, spacel2_rho_v_y, spacel2_e);
Space<double>::assign_dofs(spacesl2);
int ndofl2 = Space<double>::get_num_dofs(spacesl2);
*/


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


  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> boundary_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	



	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);


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
    EulerEquationsWeakForm_Mass wf_mass(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);


  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> * matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * dS_matrix = new CSCMatrix<double>;  
CSCMatrix<double> * dg_matrix = new CSCMatrix<double>; 
CSCMatrix<double> * K_matrix = new CSCMatrix<double>; 

 
    OGProjection<double> ogProjection;
	Lumped_Projection lumpedProjection;
/*
double* coeff_vec_init = new double[ndof];
  ogProjection.project_global(spaces,init_slns, coeff_vec_init, norms_l2 );
	Solution<double>::vector_to_solutions(coeff_vec_init, spaces, prev_slns);
Space<double>::assign_dofs(spaces);
*/

//Projection of the initial condition

	double* coeff_vec_init = new double[ndof];
  ogProjection.project_global(spaces,init_slns, coeff_vec_init, norms_l2 );
	Solution<double>::vector_to_solutions(coeff_vec_init, spaces, prev_slns);	

		KrivodonovaDiscontinuityDetector dis_detect_init(spaces, prev_slns);
		std::set<int> discont_elem =	dis_detect_init.get_discontinuous_element_ids();
		DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
		HPAdapt adapting(spaces, &error_calculator);
		bool new_adap = adapting.reduce_order(&discont_elem);

		//m1view.show(spaces[0]);
		//m2view.show(spaces[1]);m3view.show(spaces[2]);m4view.show(spaces[3]);

	 delete [] coeff_vec_init;

	//for(int i = 0; i<4;i++) spaces[i]->set_uniform_order(1);
	//Space<double>::assign_dofs(spaces);

int ndof_p1 = Space<double>::get_num_dofs(spaces);
//printf("ndof_p1=%i\n", ndof_p1);
			dof_rho = space_rho->get_num_dofs();
			dof_v_x = space_rho_v_x->get_num_dofs();
			dof_v_y = space_rho_v_y->get_num_dofs();
			dof_e = space_e->get_num_dofs();
Space<double>::assign_dofs(spaces);

			double* vec_1 = new double[ndof_p1];	double* vec_2 = new double[ndof_p1];
			double* P_plus = new double[ndof_p1]; double* P_minus = new double[ndof_p1];
			double* Q_plus = new double[ndof_p1]; double* Q_minus = new double[ndof_p1];	
			double* R_plus = new double[ndof_p1]; double* R_minus = new double[ndof_p1];

			dp_mass.set_spaces(spaces);
		  dp_mass.assemble(mass_matrix);

	bool* fct = get_fct_dofs(spaces);	
		CSCMatrix<double> * lumped_matrix = massLumping(fct,mass_matrix);


				lumpedProjection.project_lumped(spaces,init_slns, vec_1, matrix_solver, fct);
				ogProjection.project_global(spaces,init_slns, vec_2, norms_l2);
		lumped_flux_limiter(mass_matrix, lumped_matrix, vec_1, vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e,fct);


Solution<double>::vector_to_solutions(vec_1, spaces, prev_slns);

 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
delete [] vec_1; 
delete [] vec_2;


	for(int i = 0; i<4;i++) spaces[i]->set_uniform_order(P_INIT);
	Space<double>::assign_dofs(spaces);
	ndof = Space<double>::get_num_dofs(spaces);
	dof_rho = space_rho->get_num_dofs();
	dof_v_x = space_rho_v_x->get_num_dofs();
	dof_v_y = space_rho_v_y->get_num_dofs();
	dof_e = space_e->get_num_dofs();
	Space<double>::assign_dofs(spaces);
	//printf("ndof=%i\n", ndof);


	dp_mass.set_spaces(spaces);






 /*
double* coeff_vec_l2 = new double[ndofl2];
 ogProjection.project_global(spacesl2,init_slns, coeff_vec_l2, norms_l2 );
	Solution<double>::vector_to_solutions(coeff_vec_l2, spacesl2, prev_slns);*/


  MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
  MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
 MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(prev_slns));	

/*
			pressure->reinit();
			vel_x->reinit();
			vel_y->reinit();
			s1.show(prev_rho);     
			s2.show(vel_x);
			s3.show(vel_y);
  		pressure_view.show(pressure);
View::wait(HERMES_WAIT_KEYPRESS);

*/







//--------------------------------

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
	char filename[40];


	dp_mass.assemble(mass_matrix);
	mass_matrix->multiply_with_Scalar(1./time_step);


	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
 

//NumericalFlux* num_flux =new ApproxRoeNumericalFlux(KAPPA, euler_fluxes); 
//NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);
NumericalFlux* num_flux =new VijayasundaramNumericalFlux(KAPPA);


	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA,mesh, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA,mesh, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);


	EulerK wf_convection(KAPPA,  prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	EulerS wf_bdry(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,false);

	DiscreteProblem<double> dp(&wf_convection, spaces); 
	DiscreteProblem<double> dp_DG(&wf_DG, spaces);
	DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 


	EulerK wf_convection_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerS wf_bdry_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,false);

  DiscreteProblem<double> dp_init(&wf_convection_init,spaces);
  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 

		SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_bdry = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_conv = new SimpleVector<double> (ndof);
	double* coeff_vec = new double[ndof];	
	double* coeff_vec_2 = new double[ndof];

	ogProjection.project_global(spaces,prev_slns, coeff_vec, norms_l2);
	Space<double>::assign_dofs(spaces);


/*
for(int i = 0; i<ndof; i++) coeff_vec[i]= 0.;

Element* e =NULL;
		AsmList<double>*  al = new AsmList<double>;	
AsmList<double>*  al_l2 = new AsmList<double>;	
	for_all_active_elements(e, spaces[0]->get_mesh()){
		spaces[0]->get_element_assembly_list(e, al);
		spacesl2[0]->get_element_assembly_list(e, al_l2);
			for(unsigned int k = 0; k < al_l2->get_cnt(); k ++){
				 for(unsigned int j = 0; j < al->get_cnt(); j ++){
					if(al->get_idx()[j]==al_l2->get_idx()[k]){ 
						coeff_vec[al->get_dof()[j]]	= coeff_vec_l2[al_l2->get_dof()[k]];

					}
				}
			}
			


		}
Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);*/

		/*	pressure->reinit();
			vel_x->reinit();
			vel_y->reinit();
			s1.show(prev_rho);     
			s2.show(vel_x);
			s3.show(vel_y);
  		pressure_view.show(pressure);*/
			//m1view.show(spaces[0]);
//View::wait(HERMES_WAIT_KEYPRESS);
Linearizer lin_p, lin_v_x, lin_v_y, lin_rho;
		/*	lin_p.save_solution_vtk(pressure, "p_init.vtk", "pressure", false);        
			lin_v_x.save_solution_vtk(vel_x, "vx_init.vtk", "velocity_x", false); 
			lin_v_y.save_solution_vtk(vel_y, "vy_init.vtk", "velocity_y",false);     
			lin_rho.save_solution_vtk(prev_slns[0], "rho_init.vtk", "density", false);*/
int k = 1;
//Timestep loop
do
{	 
  	Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof); 	  
 	
 	 if(ts!=1){
			dp.assemble(K_matrix, vec_conv);
			dp_bdry.assemble(dS_matrix, vec_bdry);				
		  	if(DG) dp_DG.assemble(dg_matrix,vec_dg);	
		}else{
			dp_init.assemble(K_matrix, vec_conv);
			dp_bdry_init.assemble(dS_matrix, vec_bdry);
		 	if(DG) dp_DG_init.assemble(dg_matrix,vec_dg);	
		}

		//(M-theta(K+ds))u(n+1) = Sn +Ku(n) +(M-theta(Kn+ds))u(n)

		if(DG){
matrix->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
matrix->add_sparse_matrix(K_matrix);
		}else{
		matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());//L(U) = KU+SU
		}
		matrix->add_sparse_matrix(dS_matrix);
		matrix->multiply_with_Scalar(-theta);  //-theta L(U)	
		matrix->add_sparse_matrix(mass_matrix); 			//M/t - theta L(U)


	//-------------rhs: M/tau+ (1-theta)(L) u^n------------		
		vec_rhs->zero();
		if(DG) vec_rhs->add_vector(vec_dg); 
		vec_rhs->add_vector(vec_bdry); 
		vec_rhs->add_vector(vec_conv); 

	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}			

		for(int i=0; i<ndof;i++)		
					coeff_vec[i] += solver->get_sln_vector()[i];
							
		

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

	

			// Visualize the solution.


	/*	sprintf(title, "vx: ts=%i",ts);	
			s2.set_title(title);
			sprintf(title, "vy: ts=%i",ts);	
			s3.set_title(title);
			vel_x->reinit();
			vel_y->reinit();
			s2.show(vel_x);
			s3.show(vel_y);

			sprintf(title, "Pressure: ts=%i",ts);
			pressure_view.set_title(title);
			pressure->reinit();
			pressure_view.show(pressure);


			sprintf(title, "Density: ts=%i",ts);
			s1.set_title(title);
			s1.show(prev_rho);*/


	//View::wait(HERMES_WAIT_KEYPRESS);

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix->free();
		
if((current_time >= 0.1*k)&&(current_time < 0.1*k +time_step))
{
				pressure->reinit();
				vel_x->reinit();
				vel_y->reinit();
				sprintf(filename, "pressure-%i.vtk", ts );
			lin_p.save_solution_vtk(pressure, filename, "pressure", false);
        sprintf(filename, "v_x-%i.vtk", ts );
			lin_v_x.save_solution_vtk(vel_x, filename, "velocity_x", false);
  sprintf(filename, "v_y-%i.vtk", ts );
			lin_v_y.save_solution_vtk(vel_y,filename, "velocity_y",false);
       sprintf(filename, "rho-%i.vtk", ts );
			lin_rho.save_solution_vtk(prev_slns[0], filename, "density", false);
			k++;

}

}while (current_time < T_FINAL);


				pressure->reinit();
				vel_x->reinit();
				vel_y->reinit();
        
			lin_p.save_solution_vtk(pressure, "p_end.vtk", "pressure", false);        
			lin_v_x.save_solution_vtk(vel_x, "vx_end.vtk", "velocity_x", false); 
			lin_v_y.save_solution_vtk(vel_y, "vy_end.vtk", "velocity_y",false);     
			lin_rho.save_solution_vtk(prev_slns[0], "rho_end.vtk", "density", false);




		//Cleanup
			delete mass_matrix;
			delete matrix;			
delete dS_matrix;
delete dg_matrix;
delete K_matrix;




			delete[] coeff_vec_2;
			delete [] coeff_vec;
			//delete [] coeff_vec_l2;
			delete vec_rhs;
			delete vec_dg;
			delete vec_bdry;
			delete vec_conv;






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

