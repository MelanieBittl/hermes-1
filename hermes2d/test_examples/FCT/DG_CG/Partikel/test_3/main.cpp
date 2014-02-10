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

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.
const double time_step = 1e-2;
const double T_FINAL = 20;                       // Time interval length. 

const double theta = 1.;

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

  // Perform initial mesh refinements (optional).
 /* for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();

		}*/
Element* e = NULL;Node* vn=NULL;
  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();
		//y-Koord. der Knoten auf cos-Kurve verschieben
		 
			for_all_vertex_nodes(vn, basemesh)
			{	
				if(vn->bnd) 
					if((vn->x>0.)&&(vn->x<1.))		
					{		if(vn->y>0.)
							{
								
								vn->y = (Hermes::cos(vn->x*PI)+1.28/0.72)/(2.*1.28/0.72+2.) ;	
							}else if(vn->y<0.){
								vn->y =-( (Hermes::cos(vn->x*PI)+1.3/0.7)/(2.*1.3/0.7+2.) );
								//vn->y = -(Hermes::cos(vn->x*PI)+1.28/0.72)/(2.*1.28/0.72+2.) ;	
							}
					}
			}
		}






 	 mesh->copy(basemesh);

Element* test_element = RefMap::element_on_physical_coordinates(true, mesh, 1., 0);

/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;
SpaceSharedPtr<double> space_rho_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));

SpaceSharedPtr<double> space_rho_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));


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
EulerBoundary wf_bdry_low(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p, low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g, low_rho_p);

 /* EulerPenalty wf_penalty_init(SIGMA,density_particle,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
EulerPenalty wf_penalty(SIGMA, density_particle, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
*/
  EulerSource wf_source_init(density_particle, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
EulerSource wf_source(density_particle,prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
EulerSource wf_source_low(density_particle,low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g,low_rho_p);

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

DiscreteProblem<double> dp_source_low(&wf_source_low, spaces); 
DiscreteProblem<double> dp_bdry_low(&wf_bdry_low, spaces); 



//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(GAMMA);

NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(GAMMA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(GAMMA);



	EulerInterface wf_DG_init(GAMMA,mesh,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(GAMMA,mesh, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p,num_flux,euler_fluxes,riemann_invariants);


  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);



  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;
		CSCMatrix<double> matrix2;    
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
		vec_source.zero(); 

		double* coeff_vec= new double[ndof];	
		double* coeff_vec_2 = new double[ndof];




//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);




		MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach_g(new  MachNumberFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y(new VelocityFilter_y(prev_slns));
char title[100];

             /*           sprintf(title, "Pressure gas");
                        pressure_view_g.set_title(title);
                        pressure_g->reinit();
                        pressure_view_g.show(pressure_g);
                        
                        sprintf(title, "Mach gas:");
                        mach_view_g.set_title(title);
                        mach_g->reinit();
                        mach_view_g.show(mach_g);

                        sprintf(title, "density_gas");
                        s1_g.set_title(title);
                        s1_g.show(prev_rho_g);
  						sprintf(title, "v_x_gas");
                        s2_g.set_title(title);
						vel_x->reinit();
                        s2_g.show(vel_x);
  						sprintf(title, "v_y_gas");
                        s3_g.set_title(title);
						vel_y->reinit();
                        s3_g.show(vel_y);
						sprintf(title, "density_particle");
                        s1_p.set_title(title);
                        s1_p.show(prev_rho_p);*/

	//View::wait(HERMES_WAIT_KEYPRESS);
Linearizer lin;	



// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
		
double norm = 1000;
double norm_rel = 1000;
		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(&mass_matrix);
mass_matrix.multiply_with_Scalar(1./time_step);
CSCMatrix<double> * lumped_matrix;


double mach = 0.; double residual = 10000;
//Timestep loop
do
{	 
	  
 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, residual = %e", ts, current_time, ndof, residual); 		
 	// if(ts!=1){
			dp_conv.assemble(&conv_matrix, &vec_conv);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);
			dp_DG.assemble(&dg_matrix);	
				if(mach>=0.09) dp_source.assemble(&vec_source);
		
		/*}else{
			dp_conv_init.assemble(&conv_matrix, &vec_conv);
			dp_bdry_init.assemble(&bdry_matrix, &vec_bdry);
				if(mach>=0.9) dp_source.assemble(&vec_source);
			//dp_source_init.assemble(&source_matrix, &vec_source);
	
		}*/


		//matrix.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax())
matrix.create(dg_matrix.get_size(),dg_matrix.get_nnz(), dg_matrix.get_Ap(), dg_matrix.get_Ai(),dg_matrix.get_Ax());
matrix.add_sparse_matrix(&conv_matrix);
		matrix.add_sparse_matrix(&bdry_matrix);
		matrix.multiply_with_Scalar(-theta);  
		matrix.add_sparse_matrix(&mass_matrix); 



matrix2.create(bdry_matrix.get_size(),bdry_matrix.get_nnz(), bdry_matrix.get_Ap(), bdry_matrix.get_Ai(),bdry_matrix.get_Ax());	
//matrix2.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());	
		//matrix2.add_sparse_matrix(&bdry_matrix);		
		matrix2.multiply_with_Scalar(-theta); 
		matrix2.add_sparse_matrix(&mass_matrix); 
		
	//-------------rhs: ------------		
		matrix2.multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs.zero(); 
		vec_rhs.add_vector(coeff_vec_2); 		
		vec_rhs.add_vector(&vec_bdry); 
		//vec_rhs.add_vector(&vec_dg);
		//vec_rhs.add_vector(&vec_conv); 
		if(mach>=0.09) vec_rhs.add_vector(&vec_source); 

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
	coeff_vec[i] = solver->get_sln_vector()[i];

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

if(mach<0.09)
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

                        sprintf(title, "density_gas: ts=%i",ts);
                        s1_g.set_title(title);
                        s1_g.show(prev_rho_g);
  						sprintf(title, "v_x_gas: ts=%i",ts);
                        s2_g.set_title(title);
						vel_x->reinit();
                        s2_g.show(vel_x);
  						sprintf(title, "v_ye_gas: ts=%i",ts);
                        s3_g.set_title(title);
						vel_y->reinit();
                        s3_g.show(vel_y);
						sprintf(title, "density_particle: ts=%i",ts);
                        s1_p.set_title(title);
                        s1_p.show(prev_rho_p);*/



	//View::wait(HERMES_WAIT_KEYPRESS);




	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix.free();
		matrix2.free();



	if(ts%1000 ==1)
{
				pressure_g->reinit();
				mach_g->reinit();
			  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );       
			lin.save_solution_vtk(pressure_g, filename, "pressure", false);
sprintf(filename, "m-%i.vtk", ts );      
			lin.save_solution_vtk(mach_g, filename, "mach", false);

sprintf(filename, "rho_g-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[0], filename, "density_gas", false);
sprintf(filename, "rho_p-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[4], filename, "density_particle", false);

}

 



}while ((current_time < T_FINAL)&&(residual > 1e-10));




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

				pressure_g->reinit();
				mach_g->reinit();
 
			lin.save_solution_vtk(pressure_g, "pressure_end.vtk", "pressure", false);  
			lin.save_solution_vtk(mach_g, "mach_end.vtk", "mach", false); 
			lin.save_solution_vtk(prev_slns[0], "rho_g_end.vtk", "density_gas", false);  
			lin.save_solution_vtk(prev_slns[4], "rho_p_end.vtk", "density_particle", false);




			  // Clean up.

delete [] coeff_vec;
delete [] coeff_vec_2;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

