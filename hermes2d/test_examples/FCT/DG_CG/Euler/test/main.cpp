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
const int P_INIT =2;       						// Initial polynomial degree.
const double time_step = 1e-2;
const double T_FINAL = 0.5;                       // Time interval length. 

const double theta = 1.;
const bool DG = true;
// Equation parameters.  
 
     
// Kappa.
const double KAPPA = 1.4; 
    
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);
Element* e = NULL;Node* vn=NULL;

 //for_all_active_elements(e, basemesh)
	//	basemesh->refine_quad_to_triangles(e);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();

		}

 	 mesh->copy(basemesh);
double delta_x = 100; double delta_max= 0.;
for_all_active_elements(e, basemesh)
{
		for(int i = 1; i<e->get_nvert(); i++)
		{
			delta_x = std::fabs(e->vn[i]->x - e->vn[i-1]->x);
			if(delta_x>delta_max) delta_max = delta_x;
		}
}
printf("CFL = %f \n", time_step/delta_max*0.2);


/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;
/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));
/*
SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));
*/

		SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));


	int dof_rho = space_rho->get_num_dofs();
	int dof_v_x = space_rho_v_x->get_num_dofs();
	int dof_v_y = space_rho_v_y->get_num_dofs();
	int dof_e = space_e->get_num_dofs();

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new  CustomInitialCondition_rho_v_x(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new  ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));

  	MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_e(new Solution<double>);


  	MeshFunctionSharedPtr<double> diff_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_e(new Solution<double>);


  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_x(new  CustomInitialCondition_rho_v_x(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_y(new  ConstantSolution<double>(mesh, V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > diff_slns(diff_rho, diff_rho_v_x, diff_rho_v_y, diff_e);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
 ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
 ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  ScalarView mach_view("mach", new WinGeom(0, 700, 600, 300));
			s2.set_min_max_range(-1.,1.);
			s3.set_min_max_range(-1.,1.);
			mach_view.set_min_max_range(0.0, 0.1);
			//pressure_view.set_min_max_range(0.68,0.72);
			s1.set_min_max_range(0.9, 1.6);


  ScalarView s5("diff_rho", new WinGeom(700, 0, 600, 300));
ScalarView radius_view("radius_velocity", new WinGeom(0, 700, 600, 300));

radius_view.set_min_max_range(0,0.1);
Linearizer lin_p, lin_m, lin_rho, lin_r;

	/*	MeshFunctionSharedPtr<double> pressure_init(new PressureFilter(init_slns, KAPPA));


				pressure_init->reinit();
    
			lin_p.save_solution_vtk(pressure_init, "p_init.vtk", "pressure", true);


			lin_rho.save_solution_vtk(init_slns[0], "rho_init.vtk", "density", true);*/



/*	OrderView m1view("mesh", new WinGeom(1000, 0, 500, 400));m1view.show(spaces[0]);
m1view.show(spaces[0]);
View::wait(HERMES_WAIT_KEYPRESS);

/*
MeshFunctionSharedPtr<double> mach_init(new  MachNumberFilter(init_slns, KAPPA));
		s1.show(init_rho);
		mach_view.show(mach_init);
*/
//View::wait(HERMES_WAIT_KEYPRESS);

//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
 

NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);


	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

		EulerInterface wf_DG_init(KAPPA,mesh, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA,mesh, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);
	EulerK wf_convection_init(KAPPA,init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

	EulerS wf_bdry_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e,false);
	EulerS wf_bdry(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,false);





  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
  DiscreteProblem<double> dp_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp(&wf_convection, spaces); 
  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 

  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> * matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * dS_matrix = new CSCMatrix<double>;  
CSCMatrix<double> * dg_matrix = new CSCMatrix<double>; 
CSCMatrix<double> * K_matrix = new CSCMatrix<double>; 


    OGProjection<double> ogProjection;
		SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_bdry = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_conv = new SimpleVector<double> (ndof);
		double* coeff_vec = new double[ndof];	
		double* coeff_vec_2 = new double[ndof];




//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	


// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
double norm = 1000;
double norm_rel = 1000;

		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(mass_matrix);
mass_matrix->multiply_with_Scalar(1./time_step);


//Timestep loop
do
{	 
  //if(ts  % 100 == 1)
	//Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, norm =%f, norm_rel = %f", ts, current_time, ndof, norm, norm_rel); 	  
 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof); 		
 	 if(ts!=1){
			dp.assemble(K_matrix, vec_conv);
			dp_bdry.assemble(dS_matrix, vec_bdry);
		  	if(DG) dp_DG.assemble(dg_matrix,vec_dg);	
		}else{
			dp_init.assemble(K_matrix, vec_conv);
			dp_bdry_init.assemble(dS_matrix, vec_bdry);		
		 	if(DG)dp_DG_init.assemble(dg_matrix,vec_dg);	
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
		if(DG)vec_rhs->add_vector(vec_dg); 
		vec_rhs->add_vector(vec_bdry); 
		vec_rhs->add_vector(vec_conv); 


	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
		//	norm = get_l2_norm(solver->get_sln_vector(), ndof);	
		//Solution<double>::vector_to_solutions(solver->get_sln_vector(), spaces, diff_slns);	

		for(int i=0; i<ndof;i++)		
					coeff_vec[i] += solver->get_sln_vector()[i];
							


			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);


	

			// Visualize the solution.
             /*  MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
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

                        MeshFunctionSharedPtr<double> mach(new MachNumberFilter(prev_slns, KAPPA));
                        sprintf(title, "Mach: ts=%i",ts);
                        mach_view.set_title(title);
                        mach->reinit();
                        mach_view.show(mach);

 MeshFunctionSharedPtr<double> radius(new RadiusVelocityFilter(prev_slns));
radius_view.show(radius);



                        sprintf(title, "Density: ts=%i",ts);
                        s1.set_title(title);
                        s1.show(prev_rho);
						s4.show(prev_e);*/
				//s5.show(diff_slns[0]);

	//View::wait(HERMES_WAIT_KEYPRESS);

if(ts%100==1){


		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));
 MeshFunctionSharedPtr<double> radius(new RadiusVelocityFilter(prev_slns));

				pressure->reinit();
				mach->reinit();
			  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );        
			lin_p.save_solution_vtk(pressure, filename, "pressure", true);
sprintf(filename, "m-%i.vtk", ts );   
			lin_m.save_solution_vtk(mach, filename, "mach", true);
sprintf(filename, "rho-%i.vtk", ts );
			lin_rho.save_solution_vtk(prev_slns[0], filename, "density", true);
sprintf(filename, "r-%i.vtk", ts );
			lin_r.save_solution_vtk(radius, filename, "radius", true);

}
	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix->free();
 

//double abs = get_l2_norm(coeff_vec, ndof);
//norm_rel= norm/abs;

}while (current_time < T_FINAL);
//while ((current_time < T_FINAL)||(norm <1e-12)||(norm_rel<1e-08));



Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  


		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));
 MeshFunctionSharedPtr<double> radius(new RadiusVelocityFilter(prev_slns));

				pressure->reinit();
				mach->reinit();
			  char filename[40];
			  sprintf(filename, "p_end.vtk", ts );       
			lin_p.save_solution_vtk(pressure, filename, "pressure", true);
sprintf(filename, "m_end.vtk", ts );      
			lin_m.save_solution_vtk(mach, filename, "mach", true);

sprintf(filename, "rho_end.vtk", ts );  
			lin_rho.save_solution_vtk(prev_slns[0], filename, "density", true);

sprintf(filename, "r_end.vtk", ts );
			lin_r.save_solution_vtk(radius, filename, "radius", true);




		//Cleanup
			delete mass_matrix;
			delete matrix;			
delete dS_matrix;
delete dg_matrix;
delete K_matrix;



			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;
			delete vec_bdry;
			delete vec_conv;

			delete num_flux;






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

