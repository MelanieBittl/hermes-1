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
#include "fct.cpp"

const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 0.01;
const double T_FINAL = 100;                       // Time interval length. 

const double theta = 0.5;

// Equation parameters.  
 
     
// GAMMA.
const double GAMMA = 1.4; 

// Penalty Parameter.
double SIGMA = 10000000.;

//Particle density_particle
const double density_particle = 10.; 



    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;


int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain_all.mesh", basemesh);
Element* e = NULL;Node* vn=NULL;


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();

		}

 	 mesh->copy(basemesh);

Element* test_element = RefMap::element_on_physical_coordinates(true, mesh, 1, 0);

/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;

		SpaceSharedPtr<double> space_rho_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_g(new H1Space<double>(mesh, P_INIT));


		SpaceSharedPtr<double> space_rho_p(new H1Space<double>(mesh, P_INIT));	


int dof_rho = space_rho_g->get_num_dofs();
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho_g, space_rho_v_x_g, space_rho_v_y_g, space_e_g,space_rho_p);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho_g(new CustomInitialCondition_rho(mesh,  GAMMA));	
  MeshFunctionSharedPtr<double> init_rho_v_x_g(new  CustomInitialCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_v_y_g(new  ConstantSolution<double>(mesh,  0.));	
  MeshFunctionSharedPtr<double> init_rho_e_g(new CustomInitialCondition_e(mesh, GAMMA));

 MeshFunctionSharedPtr<double> init_rho_p(new ConstantSolution<double>(mesh, 0.));	


  MeshFunctionSharedPtr<double> bdry_rho_g(new BoundaryCondition_rho(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new  ConstantSolution<double>(mesh, 0.));	
  MeshFunctionSharedPtr<double> bdry_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));	

  MeshFunctionSharedPtr<double> bdry_rho_p(new ConstantSolution<double>(mesh, 0.));	



  	MeshFunctionSharedPtr<double> prev_rho_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_e_g(new Solution<double>);

  	MeshFunctionSharedPtr<double> prev_rho_p(new Solution<double>);


  	MeshFunctionSharedPtr<double> low_rho_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_x_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_y_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_e_g(new Solution<double>);

  	MeshFunctionSharedPtr<double> low_rho_p(new Solution<double>);


	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_g(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_g(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g);
	Hermes::vector<MeshFunctionSharedPtr<double> > low_slns_g(low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g);



Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
Hermes::vector<MeshFunctionSharedPtr<double> > low_slns(low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g,low_rho_p);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity

  ScalarView s1_g("rho", new WinGeom(0, 0, 400, 300));
 ScalarView s2_g("rho_v_x_g", new WinGeom(300, 0, 400, 300));
 ScalarView s3_g("rho_v_y_g", new WinGeom(300, 300, 400, 300));
ScalarView s4_g("rho_e_g", new WinGeom(300, 700, 400, 300));
  ScalarView pressure_view_g("Pressure-gas", new WinGeom(0, 300, 400, 300));
  ScalarView mach_view_g("mach-gas", new WinGeom(0, 750, 400, 300));


  ScalarView s1_p("rho_p", new WinGeom(700, 0, 400, 300));

/*
s1_g.show(init_rho_g);
s2_g.show(init_rho_v_x_g);
s3_g.show(init_rho_v_y_g);
s1_p.show(init_rho_p);


MeshFunctionSharedPtr<double> pressure_init_g(new PressureFilter(init_slns_g, GAMMA));
				pressure_view_g.show(pressure_init_g);
MeshFunctionSharedPtr<double> mach_init_g(new  MachNumberFilter(init_slns_g, GAMMA));
				mach_view_g.show(mach_init_g);*/

//----------------------Weakforms--------------------------------------------

	EulerK wf_convection_init(GAMMA,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
	EulerK wf_convection(GAMMA, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
  EulerEquationsWeakForm_Mass wf_mass;

	EulerBoundary wf_bdry_init(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
	EulerBoundary wf_bdry(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);

 /* EulerPenalty wf_penalty_init(SIGMA,density_particle,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
EulerPenalty wf_penalty(SIGMA, density_particle, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
*/
  EulerSource wf_source_init(density_particle, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
EulerSource wf_source(density_particle,prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);


  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);

  DiscreteProblem<double> dp_conv_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp_conv(&wf_convection, spaces); 

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 

/*  DiscreteProblem<double> dp_penalty_init(&wf_penalty_init,spaces); 
  DiscreteProblem<double> dp_penalty(&wf_penalty, spaces);*/ 

  DiscreteProblem<double> dp_source_init(&wf_source_init,spaces); 
  DiscreteProblem<double> dp_source(&wf_source, spaces); 



  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;
		CSCMatrix<double> matrix2;    
	CSCMatrix<double>  mass_matrix;   
	CSCMatrix<double>  source_matrix;  
	CSCMatrix<double>  conv_matrix; 
	CSCMatrix<double>  penalty_matrix;
	CSCMatrix<double>  bdry_matrix;  


    OGProjection<double> ogProjection;
	
		SimpleVector<double>  vec_rhs(ndof);
		SimpleVector<double>  vec_bdry(ndof);
		SimpleVector<double>  vec_conv(ndof);
		SimpleVector<double>  vec_source(ndof);
		SimpleVector<double>  vec_penalty(ndof);

		double* coeff_vec= new double[ndof];	
		double* coeff_vec_2 = new double[ndof];
		double* u_L= new double[ndof];;
		double* flux_vec = new double[ndof];
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	



//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);



            /*            s1_g.show(prev_rho_g);                       
                        s2_g.show(prev_rho_v_x_g);
                        s3_g.show(prev_rho_v_y_g);
                        s4_g.show(prev_rho_e_g);
                        s1_p.show(prev_rho_p);


View::wait(HERMES_WAIT_KEYPRESS);*/



		MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach_g(new  MachNumberFilter(prev_slns, GAMMA));
	
Linearizer lin;	



// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
double norm = 1000;
double norm_rel = 1000;
		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(&mass_matrix);
mass_matrix.multiply_with_Scalar(1./time_step);
CSCMatrix<double> * lumped_matrix;
lumped_matrix = massLumping(&mass_matrix);

double mach = 0.; double residual = 10000;
//Timestep loop
do
{	 
	  
 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, residual = %e", ts, current_time, ndof, residual); 		
 	 if(ts!=1){
			dp_conv.assemble(&conv_matrix, &vec_conv);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);
			//dp_source.assemble(&source_matrix, &vec_source);
				dp_source.assemble(&vec_source);
		
		}else{
			dp_conv_init.assemble(&conv_matrix, &vec_conv);
			dp_bdry_init.assemble(&bdry_matrix, &vec_bdry);
				dp_source.assemble(&vec_source);
			//dp_source_init.assemble(&source_matrix, &vec_source);
	
		}

CSCMatrix<double>* diff = NULL;	
 diff =  artificialDiffusion(GAMMA,coeff_vec,spaces,&conv_matrix);

		//matrix.create(source_matrix.get_size(),source_matrix.get_nnz(), source_matrix.get_Ap(), source_matrix.get_Ai(),source_matrix.get_Ax());
matrix.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());
//matrix.add_sparse_matrix(&dg_matrix);
		matrix.add_sparse_matrix(&bdry_matrix);
		//matrix.add_sparse_matrix(&conv_matrix);
	matrix.add_sparse_matrix(diff);

		matrix.multiply_with_Scalar(-theta); 
		matrix.add_sparse_matrix(lumped_matrix);  
		//matrix.add_sparse_matrix(&mass_matrix); 



//matrix2.create(bdry_matrix.get_size(),bdry_matrix.get_nnz(), bdry_matrix.get_Ap(), bdry_matrix.get_Ai(),bdry_matrix.get_Ax());	
matrix2.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());	
matrix2.add_sparse_matrix(diff);	
		//matrix2.add_sparse_matrix(&bdry_matrix);		
		//matrix2.multiply_with_Scalar(-theta); 
		//matrix2.add_sparse_matrix(&mass_matrix); 
		
	//-------------rhs: ------------		
		matrix2.multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs.zero(); 
		vec_rhs.add_vector(coeff_vec_2); 		
		vec_rhs.add_vector(&vec_bdry); 
		//vec_rhs.add_vector(&vec_conv); 
		if(mach>=0.1) vec_rhs.add_vector(&vec_source); 

residual = 0;
for(int i = 1; i<ndof; i++)
	residual +=vec_rhs.get(i)*vec_rhs.get(i);

	//-------------------------solution of (M-theta(K+P+B+S)) (u(n+1)-u(n) = Sn +Ku(n) +Bn+Pn------------ 
		// Hermes::Mixins::Loggable::Static::info("Solving"); 	
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (&matrix,&vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	

for(int i=0; i<ndof;i++)
	u_L[i] = solver->get_sln_vector()[i]+coeff_vec[i];



antidiffusiveFlux(&mass_matrix,lumped_matrix,diff,&matrix2, &vec_bdry,u_L, flux_vec, P_plus,P_minus, Q_plus, Q_minus,  R_plus, R_minus,dof_rho,time_step,GAMMA );
	
		for(int i=0; i<ndof;i++)		
					coeff_vec[i] = u_L[i] + flux_vec[i]*time_step/lumped_matrix->get(i,i);
							

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

if(mach<0.1)
{	double rho		= prev_rho_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0);  
	double rho_v_x	= prev_rho_v_x_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	double rho_v_y	= prev_rho_v_y_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	double rho_e	= prev_rho_e_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	mach = QuantityCalculator::calc_mach(rho, rho_v_x, rho_v_y, rho_e, GAMMA);
}

			// Visualize the solution.
 		/*Hermes::Mixins::Loggable::Static::info("Visualize"); 	
                        sprintf(title, "Pressure gas: ts=%i",ts);
                        pressure_view_g.set_title(title);
                        pressure_g->reinit();
                        pressure_view_g.show(pressure_g);
                        
                        sprintf(title, "Mach gas: ts=%i",ts);
                        mach_view_g.set_title(title);
                        mach_g->reinit();
                        mach_view_g.show(mach_g);

                        sprintf(title, "density_particle_gas: ts=%i",ts);
                        s1_g.set_title(title);
                        s1_g.show(prev_rho_g);
  						sprintf(title, "v_x_particle_gas: ts=%i",ts);
                        s2_g.set_title(title);
                        s2_g.show(prev_rho_v_x_g);
  						sprintf(title, "v_y_particle_gas: ts=%i",ts);
                        s3_g.set_title(title);
                        s3_g.show(prev_rho_v_y_g);
						sprintf(title, "density_particle_particle: ts=%i",ts);
                        s1_p.set_title(title);
                        s1_p.show(prev_rho_p);*/



	//View::wait(HERMES_WAIT_KEYPRESS);




	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix.free();
		if(diff!=NULL) delete diff;

	if(ts%100 ==1)
{
				pressure_g->reinit();
				mach_g->reinit();
			  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );       
			lin.save_solution_vtk(pressure_g, filename, "pressure", true);
sprintf(filename, "m-%i.vtk", ts );      
			lin.save_solution_vtk(mach_g, filename, "mach", true);

sprintf(filename, "rho_g-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[0], filename, "density_gas", true);
sprintf(filename, "rho_p-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[4], filename, "density_particle", true);

}

 



}while ((current_time < T_FINAL)&&(residual > 1e-10));




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

				pressure_g->reinit();
				mach_g->reinit();
 
			lin.save_solution_vtk(pressure_g, "pressure_end.vtk", "pressure", true);  
			lin.save_solution_vtk(mach_g, "mach_end.vtk", "mach", true); 
			lin.save_solution_vtk(prev_slns[0], "rho_g_end.vtk", "density_gas", true);  
			lin.save_solution_vtk(prev_slns[4], "rho_p_end.vtk", "density_particle", true);




			  // Clean up.
 			delete [] P_plus;
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;

delete [] flux_vec;
delete [] coeff_vec;
delete [] coeff_vec_2;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

