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
const double time_step = 1e-4;
const double T_FINAL = 6.;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// Kappa.
const double KAPPA = 1.4; 
// Inlet x-velocity (dimensionless).
const double V1_EXT =0.2;        
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

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();
		//y-Koord. der Knoten auf cos-Kurve verschieben
			Node* vn; 
			for_all_vertex_nodes(vn, basemesh)
			{	
				if(vn->bnd) 
					if((vn->x>0.)&&(vn->x<4.))		
					{		if(vn->y>0.)
							{
									vn->y = (Hermes::cos(PI*vn->x/2.)+3.)/4. ;
							}else if(vn->y<0.)
									vn->y = -(Hermes::cos(PI*vn->x/2.)+3.)/4. ;
					}
			}
		}
 	 mesh->copy(basemesh);




/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;
/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));*/

SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));

	/*SpaceSharedPtr<double> space_rho(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new SpaceBB<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new SpaceBB<double>(mesh, P_INIT));

		SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));*/


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
  MeshFunctionSharedPtr<double> init_rho_v_x(new  ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new  ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));

  	MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_e(new Solution<double>);


  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_x(new  ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new  ConstantSolution<double>(mesh, V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho, init_rho_v_x, init_rho_v_y, init_e);

	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
 ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
 ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
//ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  ScalarView mach_view("mach", new WinGeom(0, 700, 600, 300));
			s2.set_min_max_range(1.,4.);
			s3.set_min_max_range(-1.,1.);
			mach_view.set_min_max_range(0.17, 0.48);
			pressure_view.set_min_max_range(0.68,0.72);
			s1.set_min_max_range(0.91, 1.);

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
		NumericalFlux* num_flux = new LaxFriedrichsNumericalFlux(KAPPA);
	//NumericalFlux* num_flux = new HLLNumericalFlux(KAPPA);


	EulerInterface wf_DG_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux);
	EulerInterface wf_DG(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux);
	EulerK wf_convection_init(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass;

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
	CSCMatrix<double> * matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * mat_rhs = new CSCMatrix<double>; 
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double>;  
	CSCMatrix<double> * dS_matrix = new CSCMatrix<double>;  


    OGProjection<double> ogProjection;
		SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_in = new SimpleVector<double> (ndof);
		double* coeff_vec = new double[ndof];	
		double* coeff_vec_2 = new double[ndof];


//Projection of the initial condition
  ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	

/*
MeshFunctionSharedPtr<double> mach_2(new  MachNumberFilter(prev_slns, KAPPA));
		s1.show(prev_rho);
		mach_view.show(mach_2);*/
//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	


		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(mass_matrix);
/*
dp_bdry_init.assemble(mat_rhs, vec_in);
mat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 

for(int i =0 ; i<ndof; i++) 
if(coeff_vec_2[i]!= vec_in->get(i)) printf("%i:  %f   ,   %f \n",i, coeff_vec_2[i], vec_in->get(i));
//printf("%i, coeff = %f :  %f   ,   %f \n",i,coeff_vec[i], coeff_vec_2[i], vec_in->get(i));*/

//Timestep loop
do
{	 
  	//Hermes::Mixins::Loggable::Static::info("dofs=%i,%i,%i,%i= %i", dof_rho, dof_v_x,dof_v_y, dof_e, ndof);
  	Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof);
 	  
 	  if(ts!=1){
			//dp.assemble(mat_rhs, vec_in);
			dp.assemble(mat_rhs);
			dp_bdry.assemble(dS_matrix, vec_in);
			vec_in->multiply_with_Scalar(time_step);
		  	dp_DG.assemble(vec_dg);	
		}else{
			//dp_init.assemble(mat_rhs,vec_in);
			dp_init.assemble(mat_rhs);
			dp_bdry_init.assemble(dS_matrix, vec_in);
			vec_in->multiply_with_Scalar(time_step);
		 	dp_DG_init.assemble(vec_dg);	
		}



				matrix->create(mat_rhs->get_size(),mat_rhs->get_nnz(), mat_rhs->get_Ap(), mat_rhs->get_Ai(),mat_rhs->get_Ax());//L(U) = KU+SU
		matrix->add_sparse_matrix(dS_matrix);
		matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			matrix->add_sparse_matrix(mass_matrix); 				//M/t - theta L(U)
			mat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)|
			dS_matrix->multiply_with_Scalar(-theta*time_step); 
			mat_rhs->add_sparse_matrix(dS_matrix);
			mat_rhs->add_sparse_matrix(mass_matrix);  //M/t+(1-theta)L(U)

	//-------------rhs: M/tau+ (1-theta)(L) u^n------------		
			mat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 
			vec_dg->multiply_with_Scalar(time_step);
			vec_rhs->add_vector(vec_dg); 
			vec_rhs->add_vector(vec_in); 

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
		/*MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
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
			s1.show(prev_rho);*/


	//View::wait(HERMES_WAIT_KEYPRESS);

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix->free();

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


/*
Orderizer ord_space;
ord_space.save_orders_vtk(spaces[0], "space.vtk");*/



		//Cleanup
			delete mass_matrix;
			delete matrix;
			delete mat_rhs;

			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;

			delete num_flux;






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

