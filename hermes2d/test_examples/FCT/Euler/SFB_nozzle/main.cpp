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

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 0.4;                       // Time interval length. 

const double theta = 0.5;

// Equation parameters.  
 
  const double R = 287;    
// GAMMA.
const double GAMMA = 1.4;  
bool view_3D = false;
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
  mloader.load("domain_all_inlets.mesh", basemesh);

  /*MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(basemesh);
   View::wait();*/
  
  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);
   
  /*MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();
*/

	SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));	

		int dof_rho = space_rho->get_num_dofs();
		int dof_v_x = space_rho_v_x->get_num_dofs();
		int dof_v_y = space_rho_v_y->get_num_dofs();
		int dof_e = space_e->get_num_dofs();

    Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
		Space<double>::assign_dofs(spaces);
  	int ndof = Space<double>::get_num_dofs(spaces);
  	printf("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
		MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh,  GAMMA));	
		MeshFunctionSharedPtr<double> init_rho_v_x(new   CustomInitialCondition_rho_v_x(mesh, GAMMA));	
		MeshFunctionSharedPtr<double> init_rho_v_y(new   CustomInitialCondition_rho_v_y(mesh, GAMMA));	
		MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,GAMMA));	

  MeshFunctionSharedPtr<double> bdry_rho_g(new BoundaryCondition_rho(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new  BoundaryCondition_rho_v_y(mesh, GAMMA));	
  MeshFunctionSharedPtr<double> bdry_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));	



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
  ScalarView s1_g("rho", new WinGeom(0, 0, 400, 300));
 ScalarView s2_g("rho_v_x_g", new WinGeom(300, 0, 400, 300));
 ScalarView s3_g("rho_v_y_g", new WinGeom(300, 300, 400, 300));
ScalarView s4_g("rho_e_g", new WinGeom(300, 700, 400, 300));
  ScalarView pressure_view_g("Pressure-gas", new WinGeom(0, 300, 400, 300));
  ScalarView mach_view_g("mach-gas", new WinGeom(0, 750, 400, 300));
  ScalarView temp_view_g("Temperature-gas", new WinGeom(200, 750, 400, 300));

//--------------Weakforms------------
  EulerEquationsWeakForm_K  wf_K_init(GAMMA, init_slns);
  EulerBoundary wf_boundary_init(GAMMA, Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e,bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g),false );
  EulerEquationsWeakForm_Mass wf_mass;

  EulerEquationsWeakForm_K*  wf_K= new EulerEquationsWeakForm_K(GAMMA, prev_slns);
  EulerBoundary* wf_boundary= new EulerBoundary(GAMMA,
Hermes::vector<MeshFunctionSharedPtr<double> > (prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g),false);
  EulerEquationsWeakForm_K*  wf_K_low = new EulerEquationsWeakForm_K(GAMMA, low_slns);
  EulerBoundary* wf_boundary_low = new EulerBoundary(GAMMA,Hermes::vector<MeshFunctionSharedPtr<double> > (low_rho, low_rho_v_x, low_rho_v_y, low_e, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g ),false);


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
		SimpleVector<double> * vec_bdry = new SimpleVector<double>(ndof);
	double* u_L = NULL; 

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	

	double* dt_u_L = new double[ndof];

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
		MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach_g(new  MachNumberFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y(new VelocityFilter_y(prev_slns));
		MeshFunctionSharedPtr<double> temp_g(new TempFilter(prev_slns,R, GAMMA));
Linearizer lin;	

				pressure_g->reinit();
				mach_g->reinit();		
				temp_g->reinit();
				//temp_view_g.show(temp_g);


//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];

mass_matrix->multiply_with_Scalar(1./time_step);
lumped_matrix->multiply_with_Scalar(1./time_step);

Space<double>::assign_dofs(spaces);
//Timestep loop
do
{	 
  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time); 	  
 	  if(ts!=1) {
			dp_boundary->assemble(matrix_dS,vec_bdry);
		  dp_K->assemble(lowmat_rhs);
		}else{
			dp_boundary_init.assemble(matrix_dS,vec_bdry);
		  dp_K_init.assemble(lowmat_rhs);				
		}
	//------------------------artificial DIFFUSION D---------------------------------------	
			CSCMatrix<double> * diffusion = artificialDiffusion(GAMMA,coeff_vec,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);

			lowmat_rhs->add_sparse_matrix(diffusion); //L(U)=K+D
			//lowmat_rhs->add_sparse_matrix(matrix_dS); //L(U)+dS(U) 

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->add_sparse_matrix(matrix_dS); 
			low_matrix->multiply_with_Scalar(-theta);  //-theta L(U)
			low_matrix->add_sparse_matrix(lumped_matrix); 				//M_L/t - theta L(U)

			//lowmat_rhs->multiply_with_Scalar((1.0-theta));  //(1-theta)L(U)
			//lowmat_rhs->add_sparse_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)

	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 
			low_matrix->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->add_vector(coeff_vec_2); 
				vec_rhs->add_vector(vec_bdry); 
				
	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
			try{
			 lowOrd->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}		
			u_L = lowOrd->get_sln_vector();  
      Solution<double>::vector_to_solutions(u_L, spaces,low_slns);

			dp_boundary_low->assemble(vec_bdry);	
    	dp_K_low->assemble(matrix_L_low);
			CSCMatrix<double> * diffusion_low = artificialDiffusion(GAMMA,u_L,spaces, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_sparse_matrix(diffusion_low); //L(U)
			//matrix_L_low->add_sparse_matrix(matrix_dS_low); //L(U)+dS(U) 
			matrix_L_low->multiply_with_vector(u_L, dt_u_L);
			for(int i=0; i<ndof;i++)
							dt_u_L[i] +=vec_bdry->get(i);

		//---------------------------------------antidiffusive fluxes-----------------------------------	
		//Hermes::Mixins::Loggable::Static::info("antidiffusive fluxes ");
		antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, dt_u_L, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

		for(int i=0; i<ndof;i++)
					coeff_vec[i] = u_L[i]+ coeff_vec_2[i]/lumped_matrix->get(i,i);					

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);	
	

	// Visualize the solution.
 	/*	Hermes::Mixins::Loggable::Static::info("Visualize"); 	
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
                        s1_g.show(prev_rho);
  						sprintf(title, "v_x_gas: ts=%i",ts);
                        s2_g.set_title(title);
						vel_x->reinit();
                        s2_g.show(vel_x);
  						sprintf(title, "v_ye_gas: ts=%i",ts);
                        s3_g.set_title(title);
						vel_y->reinit();
                        s3_g.show(vel_y);
					
								
								temp_g->reinit();
								sprintf(title, "Temp gas: ts=%i",ts);
                        temp_view_g.set_title(title);
								 temp_view_g.show(temp_g);

*/
  		

		//	View::wait(HERMES_WAIT_KEYPRESS);



// Visualization.
    if((ts - 1) % 1000 == 0) 
    {

        pressure_g->reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-%i.vtk", ts - 1);
        lin_pressure.save_solution_vtk(pressure_g, filename, "Pressure", view_3D);
		   mach_g->reinit();
		  			        Linearizer lin_mach_t;
							  sprintf(filename, "mach-%i.vtk", ts - 1);
			lin_mach_t.save_solution_vtk(mach_g, filename, "mach", view_3D);
					   temp_g->reinit();
		  			        Linearizer lin_t;
							  sprintf(filename, "temp-%i.vtk", ts - 1);
			lin_t.save_solution_vtk(temp_g, filename, "temp", view_3D);
			
			        Linearizer lin_v_x_t;
					   sprintf(filename, "vx-%i.vtk", ts - 1);
			lin_v_x_t.save_solution_vtk(vel_x, filename, "velocity_x", view_3D);
        Linearizer lin_v_y_t;
		    sprintf(filename, "vy-%i.vtk", ts - 1);
			lin_v_y_t.save_solution_vtk(vel_y, filename, "velocity_y",view_3D);
        Linearizer lin_rho_t;
		    sprintf(filename, "rho-%i.vtk", ts - 1);
			lin_rho_t.save_solution_vtk(prev_slns[0], filename, "density", view_3D);

      
    }




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

				pressure_g->reinit();
				vel_x->reinit();
				vel_y->reinit();
 mach_g->reinit();
temp_g->reinit();
        Linearizer lin_p;
			lin_p.save_solution_vtk(pressure_g, "p_end.vtk", "pressure", view_3D);
			        Linearizer lin_mach;
			lin_mach.save_solution_vtk(mach_g, "mach_end.vtk", "mach", view_3D);
        Linearizer lin_temp;
			lin_temp.save_solution_vtk(temp_g, "temp_end.vtk", "temp",view_3D);
        Linearizer lin_v_x;
			lin_v_x.save_solution_vtk(vel_x, "vx_end.vtk", "velocity_x", view_3D);
        Linearizer lin_v_y;
			lin_v_y.save_solution_vtk(vel_y, "vy_end.vtk", "velocity_y",view_3D);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_slns[0], "rho_end.vtk", "density", view_3D);
			


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

