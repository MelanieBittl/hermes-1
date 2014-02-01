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

#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 1e-5;
const double T_FINAL = 50000000000;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// GAMMA.
const double GAMMA = 1.4; 

// Penalty Parameter.
const double SIGMA = std::pow(10,8.); 

//Particle density
const double density = 4000.; 
    

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


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();

		}

 	 mesh->copy(basemesh);


/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;
/*
		SpaceSharedPtr<double> space_rho_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_g(new H1Space<double>(mesh, P_INIT));


		SpaceSharedPtr<double> space_rho_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_p(new H1Space<double>(mesh, P_INIT));
*/

SpaceSharedPtr<double> space_rho_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));

SpaceSharedPtr<double> space_rho_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));


//int dof_1 = space_rho_g->get_num_dofs();

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho_g, space_rho_v_x_g, space_rho_v_y_g, space_e_g,space_rho_p, space_rho_v_x_p, space_rho_v_y_p, space_e_p);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho_g(new ConstantSolution<double>(mesh,  6.));	
  MeshFunctionSharedPtr<double> init_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_v_y_g(new  ConstantSolution<double>(mesh,  0.));	
  MeshFunctionSharedPtr<double> init_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));

 MeshFunctionSharedPtr<double> init_rho_p(new ConstantSolution<double>(mesh, 3.4));	
  MeshFunctionSharedPtr<double> init_rho_v_x_p(new  BoundaryCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_v_y_p(new  ConstantSolution<double>(mesh,  0.));	
  MeshFunctionSharedPtr<double> init_rho_e_p(new BoundaryCondition_rho_e(mesh, GAMMA,true));


  MeshFunctionSharedPtr<double> bdry_rho_g(new ConstantSolution<double>(mesh,  6.));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new  ConstantSolution<double>(mesh, 0.));	
  MeshFunctionSharedPtr<double> bdry_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));	

  MeshFunctionSharedPtr<double> bdry_rho_p(new ConstantSolution<double>(mesh, 3.4));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_p(new  BoundaryCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_p(new  ConstantSolution<double>(mesh, 0.));	
  MeshFunctionSharedPtr<double> bdry_rho_e_p(new BoundaryCondition_rho_e(mesh, GAMMA,true));	


  	MeshFunctionSharedPtr<double> prev_rho_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_e_g(new Solution<double>);

  	MeshFunctionSharedPtr<double> prev_rho_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_e_p(new Solution<double>);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_g(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_g(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_p(prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_p(init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);


Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);


	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity

  ScalarView s1_g("rho", new WinGeom(0, 0, 400, 300));
 ScalarView s2_g("rho_v_x_g", new WinGeom(300, 0, 400, 300));
 ScalarView s3_g("rho_v_y_g", new WinGeom(300, 300, 400, 300));
ScalarView s4_g("rho_e_g", new WinGeom(300, 700, 400, 300));
  ScalarView pressure_view_g("Pressure-gas", new WinGeom(0, 300, 400, 300));
  ScalarView mach_view_g("mach-gas", new WinGeom(0, 750, 400, 300));


  ScalarView s1_p("rho_p", new WinGeom(700, 0, 400, 300));
 ScalarView s2_p("rho_v_x_p", new WinGeom(900, 0, 400, 300));
 ScalarView s3_p("rho_v_y_p", new WinGeom(900, 300, 400, 300));
ScalarView s4_p("rho_e_p", new WinGeom(900, 700, 400, 300));

  ScalarView pressure_view_p("Pressure_p", new WinGeom(700, 300, 400, 300));
  ScalarView mach_view_p("mach_p", new WinGeom(700, 750, 400, 300));



s1_g.show(init_rho_g);
s1_p.show(init_rho_p);
s2_p.show(init_rho_v_x_p);
s3_p.show(init_rho_v_y_p);
s4_p.show(init_rho_e_p);


MeshFunctionSharedPtr<double> pressure_init_g(new PressureFilter(init_slns_g, GAMMA));
				pressure_view_g.show(pressure_init_g);
MeshFunctionSharedPtr<double> mach_init_g(new  MachNumberFilter(init_slns_g, GAMMA));
				mach_view_g.show(mach_init_g);





//View::wait(HERMES_WAIT_KEYPRESS);

//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(GAMMA);

NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(GAMMA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(GAMMA);


/*
	EulerInterface wf_DG_init(GAMMA,mesh,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(GAMMA,mesh, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p,num_flux,euler_fluxes,riemann_invariants);


  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);

*/

	EulerK wf_convection_init(GAMMA,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
	EulerK wf_convection(GAMMA, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
  EulerEquationsWeakForm_Mass wf_mass;

	EulerBoundary wf_bdry_init(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
	EulerBoundary wf_bdry(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

  EulerPenalty wf_penalty_init(SIGMA,density,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
EulerPenalty wf_penalty(SIGMA, density, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

  EulerSource wf_source_init(density, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
EulerSource wf_source(density, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);


  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);

  DiscreteProblem<double> dp_conv_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp_conv(&wf_convection, spaces); 

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 

  DiscreteProblem<double> dp_penalty_init(&wf_penalty_init,spaces); 
  DiscreteProblem<double> dp_penalty(&wf_penalty, spaces); 

  DiscreteProblem<double> dp_source_init(&wf_source_init,spaces); 
  DiscreteProblem<double> dp_source(&wf_source, spaces); 



  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;  
	CSCMatrix<double>  mass_matrix;   
	CSCMatrix<double>  source_matrix;  
	CSCMatrix<double>  conv_matrix; 
	CSCMatrix<double>  penalty_matrix;
	CSCMatrix<double>  bdry_matrix;  
	CSCMatrix<double> dg_matrix; 



    OGProjection<double> ogProjection;
		SimpleVector<double>  vec_dg(ndof);
		SimpleVector<double>  vec_rhs(ndof);
		SimpleVector<double>  vec_bdry(ndof);
		SimpleVector<double>  vec_conv(ndof);
		SimpleVector<double>  vec_source(ndof);
		SimpleVector<double>  vec_penalty(ndof);

		double coeff_vec[ndof];	
		double coeff_vec_2[ndof];




//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);



            /*            s1_g.show(prev_rho_g);                       
                        s2_g.show(prev_rho_v_x_g);
                        s3_g.show(prev_rho_v_y_g);
                        s4_g.show(prev_rho_e_g);
                        s1_p.show(prev_rho_p);
                        s2_p.show(prev_rho_v_x_p);
                        s3_p.show(prev_rho_v_y_p);
                        s4_p.show(prev_rho_e_p);


View::wait(HERMES_WAIT_KEYPRESS);*/

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
double norm = 1000;
double norm_rel = 1000;
 Hermes::Mixins::Loggable::Static::info("Assembling mass"); 	
		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(&mass_matrix);
mass_matrix.multiply_with_Scalar(1./time_step);
CSCMatrix<double> * lumped_matrix;
//lumped_matrix = massLumping(&mass_matrix);
//Timestep loop
do
{	 
	  
 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof); 		
 	 if(ts!=1){
			dp_conv.assemble(&conv_matrix, &vec_conv);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);
			dp_penalty.assemble(&penalty_matrix, &vec_penalty);
			dp_source.assemble(&source_matrix, &vec_source);
				
		  	//dp_DG.assemble(&dg_matrix,&vec_dg);	
		}else{
			dp_conv_init.assemble(&conv_matrix, &vec_conv);
			dp_bdry_init.assemble(&bdry_matrix, &vec_bdry);
			dp_penalty_init.assemble(&penalty_matrix, &vec_penalty);
			dp_source_init.assemble(&source_matrix, &vec_source);
			
		 	//dp_DG_init.assemble(&dg_matrix,&vec_dg);	
		}

		

//matrix.create(dg_matrix.get_size(),dg_matrix.get_nnz(), dg_matrix.get_Ap(), dg_matrix.get_Ai(),dg_matrix.get_Ax());
//matrix.add_sparse_matrix(&source_matrix);
//CSCMatrix<double>* diff =  artificialDiffusion(GAMMA,coeff_vec,spaces,&conv_matrix);

		matrix.create(source_matrix.get_size(),source_matrix.get_nnz(), source_matrix.get_Ap(), source_matrix.get_Ai(),source_matrix.get_Ax());
//matrix.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());
//matrix.add_sparse_matrix(&dg_matrix);
		matrix.add_sparse_matrix(&bdry_matrix);
		matrix.add_sparse_matrix(&conv_matrix);
//matrix.add_sparse_matrix(diff);
		matrix.add_sparse_matrix(&penalty_matrix);
		matrix.multiply_with_Scalar(-theta); 
//matrix.add_sparse_matrix(lumped_matrix);  
		matrix.add_sparse_matrix(&mass_matrix); 			

	//-------------rhs: ------------		
	//	matrix.multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs.zero(); 
		//vec_rhs.add_vector(coeff_vec_2); 
		//vec_rhs.add_vector(&vec_dg); 
		vec_rhs.add_vector(&vec_bdry); 
		vec_rhs.add_vector(&vec_conv); 
		vec_rhs.add_vector(&vec_penalty); 
		vec_rhs.add_vector(&vec_source); 


	//-------------------------solution of (M-theta(K+P+B+S)) (u(n+1)-u(n) = Sn +Ku(n) +Bn+Pn------------ 
 Hermes::Mixins::Loggable::Static::info("Solving"); 		
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (&matrix,&vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
 Hermes::Mixins::Loggable::Static::info("Solved"); 		
		for(int i=0; i<ndof;i++)		
					coeff_vec[i] += solver->get_sln_vector()[i];
							

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);


	

			// Visualize the solution.
 Hermes::Mixins::Loggable::Static::info("Visualize"); 	
                        MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns_g, GAMMA));
                        sprintf(title, "Pressure gas: ts=%i",ts);
                        pressure_view_g.set_title(title);
                        pressure_g->reinit();
                        pressure_view_g.show(pressure_g);

                        MeshFunctionSharedPtr<double> mach_g(new MachNumberFilter(prev_slns_g, GAMMA));
                        sprintf(title, "Mach gas: ts=%i",ts);
                        mach_view_g.set_title(title);
                        mach_g->reinit();
                        mach_view_g.show(mach_g);


                        sprintf(title, "Density_gas: ts=%i",ts);
                        s1_g.set_title(title);
                        s1_g.show(prev_rho_g);

						sprintf(title, "Density_particle: ts=%i",ts);
                        s1_p.set_title(title);
                        s1_p.show(prev_rho_p);

sprintf(title, "vel_x_particle: ts=%i",ts);
                        s2_p.set_title(title);
                        s2_p.show(prev_rho_v_x_p);

sprintf(title, "vel_y_particle: ts=%i",ts);
                        s3_p.set_title(title);
                        s3_p.show(prev_rho_v_y_p);

sprintf(title, "energy_particle: ts=%i",ts);
                        s4_p.set_title(title);
                        s4_p.show(prev_rho_e_p);


	//View::wait(HERMES_WAIT_KEYPRESS);




	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix.free();
		//delete diff;
 



}while (current_time < T_FINAL);




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

/*
		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, GAMMA));
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

*/







  // Wait for the view to be closed.
  View::wait();
  return 0;
}

