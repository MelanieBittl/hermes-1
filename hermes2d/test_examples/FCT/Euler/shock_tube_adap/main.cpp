#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "prev_solution.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;



const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 0.231;                       // Time interval length. 
const double EPS = 1e-8;
const double theta = 0.5;

const double h_max = 0.1;
const int NDOF_STOP = 6000000;
 
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;  

const double THRESHOLD = 0.3;

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 



//FCT & p-Adaptivity
#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"
#include "h_adapt.cpp"
#include "z_z.cpp"

//Visualization
const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION =false;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 500;													//Every VTK_FREQth time step the solution is saved as VTK output. 
     
const int UNREF_FREQ = 5;                         // Every UNREF_FREQth time step the mesh is derefined.
const int  UNREF_METHOD =1; 



int main(int argc, char* argv[])
{
   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();
 	 mesh->copy(basemesh);

	 double diag;			
	Element* e =NULL;
	for_all_active_elements(e, mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8.;

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

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(mesh,KAPPA));	

    MeshFunctionSharedPtr<double> prev_rho(new PrevSolution);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new PrevSolution);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new PrevSolution);
    MeshFunctionSharedPtr<double> prev_e(new PrevSolution);

    MeshFunctionSharedPtr<double> new_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> new_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> new_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> new_e(new Solution<double>);

	
  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(mesh));	
  MeshFunctionSharedPtr<double> boundary_v_x(new   ConstantSolution<double>(mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new   ConstantSolution<double>(mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(mesh,KAPPA));	
	
    MeshFunctionSharedPtr<double> low_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double>	low_rho_e(new Solution<double>);

  MeshFunctionSharedPtr<double>  R_h_1(new Solution<double>);
  MeshFunctionSharedPtr<double>  R_h_2(new Solution<double>);

  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
Hermes::vector<MeshFunctionSharedPtr<double> > low_slns(low_rho,low_rho_v_x,low_rho_v_y,low_rho_e);
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

//------------


  EulerEquationsWeakForm_K  wf_K_init(KAPPA, init_slns);
  EulerBoundary wf_boundary_init(KAPPA, Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e) );

  EulerEquationsWeakForm_Mass wf_mass;
  EulerEquationsWeakForm_K  wf_K(KAPPA, prev_slns);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, low_slns);
  EulerBoundary wf_boundary(KAPPA,
Hermes::vector<MeshFunctionSharedPtr<double> > (prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e));
 EulerBoundary wf_boundary_low(KAPPA, Hermes::vector<MeshFunctionSharedPtr<double> > (low_rho, low_rho_v_x, low_rho_v_y, low_rho_e, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e ));



  // Set up the solver, matrix, and rhs according to the solver selection.
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;   
	CSCMatrix<double> * matrix_L_low = new CSCMatrix<double> ; 
	CSCMatrix<double> * low_matrix = new CSCMatrix<double> ;  
	CSCMatrix<double> * lowmat_rhs = new CSCMatrix<double> ; 
	CSCMatrix<double> * matrix_dS = new CSCMatrix<double> ; 
	CSCMatrix<double> * matrix_dS_low = new CSCMatrix<double> ; 

				double* u_L = NULL; 
    OGProjection<double> ogProjection;
		Lumped_Projection lumpedProjection;

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;	int as = 1;
		char title[100]; char filename[40];
	
	AsmList<double>*  al = new AsmList<double>;	
	 DiscreteProblem<double> *dp_boundary; 
   DiscreteProblem<double>* dp_K;

  DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);
  AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);

//Timestep loop
do
{	 
 	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time);


	if ((ts > 1 && ts % UNREF_FREQ == 0)||(space_rho->get_num_dofs() >= NDOF_STOP)) 
    { 
    	Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1:  	 mesh->copy(basemesh);
								space_rho->set_mesh(mesh);			space_rho->set_uniform_order(P_INIT); 
								space_rho_v_x->set_mesh(mesh);	space_rho_v_x->set_uniform_order(P_INIT);  
								space_rho_v_y->set_mesh(mesh);	space_rho_v_y->set_uniform_order(P_INIT);  
								space_e->set_mesh(mesh);	space_e->set_uniform_order(P_INIT);   
                break;
        case 2: mesh->unrefine_all_elements();
								space_rho->set_mesh(mesh);			space_rho->set_uniform_order(P_INIT); 
								space_rho_v_x->set_mesh(mesh);	space_rho_v_x->set_uniform_order(P_INIT);  
								space_rho_v_y->set_mesh(mesh);	space_rho_v_y->set_uniform_order(P_INIT);  
								space_e->set_mesh(mesh);	space_e->set_uniform_order(P_INIT);   
                break;
        case 3: mesh->unrefine_all_elements();
								space_rho->set_mesh(mesh);			space_rho->set_uniform_order(P_INIT); 
								space_rho_v_x->set_mesh(mesh);	space_rho_v_x->set_uniform_order(P_INIT);  
								space_rho_v_y->set_mesh(mesh);	space_rho_v_y->set_uniform_order(P_INIT);  
								space_e->set_mesh(mesh);	space_e->set_uniform_order(P_INIT);   
                break;
        default: Exceptions::Exception("Wrong global derefinement method.");
      }
       	       
    }
    
		as=1;

	//Adaptivity loop
	do{
	  	Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f, adap: %i", ts, current_time,as);

	Space<double>::assign_dofs(spaces);

	ndof = Space<double>::get_num_dofs(spaces);

		dof_rho = space_rho->get_num_dofs();
		dof_v_x = space_rho_v_x->get_num_dofs();
		dof_v_y = space_rho_v_y->get_num_dofs();
		dof_e = space_e->get_num_dofs();
	

			double* coeff_vec = new double[ndof];
			double* coeff_vec_2 = new double[ndof];
			SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);

			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
			double* Q_plus_old = new double[ndof]; double* Q_minus_old = new double[ndof];	
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	


  // Initialize the FE problem.
    DiscreteProblem<double>* dp_mass = new DiscreteProblem<double>(&wf_mass, spaces);
  if(ts==1){
 	 dp_boundary = new DiscreteProblem<double>(&wf_boundary_init, spaces);
 		dp_K = new DiscreteProblem<double>(&wf_K_init, spaces);	  
  }else{
 		dp_boundary = new DiscreteProblem<double>(&wf_boundary, spaces);
 		dp_K = new DiscreteProblem<double>(&wf_K, spaces);		  
		}

	 Space<double>::assign_dofs(spaces); 

    dp_mass->assemble(mass_matrix);

	//----------------------MassLumping M_L--------------------------------------------------------------------
						// 	Hermes::Mixins::Loggable::Static::info("mass-lumping");
		CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of previous timestep solution / initial data

if(ts==1){
	lumpedProjection.project_lumped(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec, matrix_solver);
   ogProjection.project_global(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(init_rho, init_rho_v_x, init_rho_v_y, init_e), coeff_vec_2,  Hermes::vector<NormType> (HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM));
   }else{
   	lumpedProjection.project_lumped(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), coeff_vec, matrix_solver);
   ogProjection.project_global(spaces,Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), coeff_vec_2, Hermes::vector<NormType> (HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM));   
   }

		lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e);

Space<double>::assign_dofs(spaces);

			dp_boundary->assemble(matrix_dS);
		  dp_K->assemble(lowmat_rhs);

					//------------------------artificial DIFFUSION D---------------------------------------		
					// 	Hermes::Mixins::Loggable::Static::info("artificial Diffusion");
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
						//  	Hermes::Mixins::Loggable::Static::info("Solution low order ");
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);	
  try
  {
   lowOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(low_rho,low_rho_v_x,low_rho_v_y,low_rho_e));

	
	Space<double>::assign_dofs(spaces);
	
	if(as==1){
								//Hermes::Mixins::Loggable::Static::info("Adapt-step 1");
		HPAdapt * adapting = new HPAdapt(spaces, &error_calculator);
		int* elements_to_refine = new int[space_rho->get_mesh()->get_max_element_id()]; 
		int* no_of_refinement_steps = new int[space_rho->get_mesh()->get_max_element_id()];	
		double* elem_error = new double[space_rho->get_mesh()->get_max_element_id()];
						 // 	Hermes::Mixins::Loggable::Static::info("calc error");
		for_all_active_elements(e, space_rho->get_mesh()){	elements_to_refine[e->id] = 2; no_of_refinement_steps[e->id]=0;	}
		calc_elem_error(space_rho, low_rho,R_h_1, R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error);
		//calc_elem_error(space_rho_v_x, low_rho_v_x,R_h_1, R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error);
		//calc_elem_error(space_rho_v_y, low_rho_v_y,R_h_1, R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error);
		//calc_elem_error(space_e, low_rho_e,R_h_1, R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error);
		adapting->adapt(elements_to_refine,no_of_refinement_steps, h_min,h_max,NDOF_STOP);
	Space<double>::assign_dofs(spaces);
			delete adapting;
			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
			delete [] elem_error;			

	}else{

 //	Hermes::Mixins::Loggable::Static::info("Adapt-step 2");
			 	 DiscreteProblem<double>* dp_boundary_low = new DiscreteProblem<double>(&wf_boundary_low, spaces);
				DiscreteProblem<double>* dp_K_low = new DiscreteProblem<double>(&wf_K_low, spaces);	

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

					Solution<double>::vector_to_solutions(coeff_vec, spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(new_rho, new_rho_v_x, new_rho_v_y, new_e));	

				 prev_rho->copy(new_rho); 
				 prev_rho->set_own_mesh(new_rho->get_mesh());
				 prev_rho_v_x->copy(new_rho_v_x); 
				 prev_rho_v_x->set_own_mesh(new_rho_v_x->get_mesh());
				 prev_rho_v_y->copy(new_rho_v_y); 
				 prev_rho_v_y->set_own_mesh(new_rho_v_y->get_mesh());
				 prev_e->copy(new_e); 
				 prev_e->set_own_mesh(new_e->get_mesh()); 
				 

				if(HERMES_VISUALIZATION){
		PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
	VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
	  VelocityFilter vel_y(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 2);

					sprintf(title, "pressure: ts=%i",ts);
					pressure_view.set_title(title);
					s1.show(prev_rho);
				s2.show(&vel_x);
					s3.show(&vel_y);
					pressure_view.show(&pressure);

	  		}

			  		      // Output solution in VTK format.
/*		 if((VTK_VISUALIZATION)&&(ts  % VTK_FREQ == 0)) 
		 {
		 		PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
	VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
    	Linearizer lin_p,lin_rho,lin_v_x; 	Orderizer ord_p;		
      sprintf(filename, "pressure-%i.vtk", ts);
      lin_p.save_solution_vtk(&pressure, filename, "Pressure", true);  
      sprintf(filename, "density-%i.vtk", ts);
      lin_rho.save_solution_vtk(prev_rho, filename, "density", true); 
      sprintf(filename, "mesh-%i.vtk", ts); 
     	ord_p.save_mesh_vtk(space_rho, filename);  
         sprintf(filename, "vel_x-%i.vtk", ts); 	
     		lin_v_x.save_solution_vtk(&vel_x,filename, "velocity_x", true);

		 }*/

					delete diffusion_low;
					delete dp_boundary_low;
					delete dp_K_low;	
	}


			  // Clean up.
		delete lowOrd;
		delete diffusion;
		low_matrix->free();
		delete lumped_matrix;
		delete vec_rhs;
		delete [] P_plus;
		delete [] P_minus;
		delete [] Q_plus;
		delete [] Q_minus;
		delete [] R_plus;
		delete [] R_minus;
		delete [] coeff_vec_2;
		delete [] coeff_vec; 
		delete dp_mass;
		delete dp_boundary;
		delete dp_K;


		as++;

	}while(as<3);


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



}
while (current_time < T_FINAL);


 	/*if(VTK_VISUALIZATION)
     {
       PressureFilter pressure(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 2);
        Linearizer lin_p, lin_v_x, lin_v_y, lin_rho; 	Orderizer ord_p;		
			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
			
		Linearizer lin_rho; 	Orderizer ord_p;	
			lin_rho.save_solution_vtk(prev_rho, "rho_end.vtk", "density", true);
			     ord_p.save_mesh_vtk(space_rho, "mesh_end.vtk"); 
		}
 */

		//Cleanup
			delete mass_matrix;
			delete matrix_L_low;
			delete matrix_dS;			
			delete matrix_dS_low;	
			delete low_matrix;
			delete lowmat_rhs;



  // Wait for the view to be closed.
  View::wait();
  return 0;
}

