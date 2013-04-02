#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "boundary_condition.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;



const int INIT_REF_NUM= 5;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX =1;
const double h_max = 0.1;
const int NDOF_STOP = 100000;
const double EPS = 1e-8;
const double time_step =2e-3;
const double T_FINAL = 0.13;                       // Time interval length. 

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
#include "h_adapt.cpp"
#include "z_z.cpp"


//Visualization
const bool HERMES_VISUALIZATION = false;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool HERMES_VISUALIZATION_LOW = false;           // Set to "false" to suppress Hermes OpenGL visualization for low_order.
const bool HERMES_VISUALIZATION_ORDER = false;           // Set to "false" to suppress Hermes OpenGL visualization for low_order.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 500;													//Every VTK_FREQth time step the solution is saved as VTK output. 
     
const int UNREF_FREQ = 50;                         // Every UNREF_FREQth time step the mesh is derefined.
const int  UNREF_METHOD =1;

const bool FCT_W_HANGING_NODES = true; 

int main(int argc, char* argv[])
{
   // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

	 double diag;			
	Element* e =NULL;
	for_all_active_elements(e, &mesh){diag = e->get_diameter(); break;}
	const double h_min = diag/8.;
	info("diag=%f, h_min = %f", diag, h_min);

	H1Space<double> space_rho(&mesh,P_INIT);
  H1Space<double> space_rho_v_x(&mesh,P_INIT);
  H1Space<double> space_rho_v_y(&mesh,P_INIT);
  H1Space<double> space_e(&mesh, P_INIT);

	int dof_rho = space_rho.get_num_dofs();
	int dof_v_x = space_rho_v_x.get_num_dofs();
	int dof_v_y = space_rho_v_y.get_num_dofs();
	int dof_e = space_e.get_num_dofs();


  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
	CustomInitialCondition_rho prev_rho(&mesh);
  ConstantSolution<double> prev_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e prev_e(&mesh, KAPPA);

 Solution<double>  new_rho,new_rho_v_x,new_rho_v_y,new_e;
 Solution<double> high_rho,high_rho_v_x,high_rho_v_y,high_rho_e;

//InitialCondition auch fuer low setzen wegen Filtern!
	CustomInitialCondition_rho low_rho(&mesh);
  ConstantSolution<double> low_rho_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> low_rho_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e low_rho_e(&mesh, KAPPA);

	CustomInitialCondition_rho boundary_rho(&mesh);
  ConstantSolution<double> boundary_v_x(&mesh,  V1_EXT);
  ConstantSolution<double> boundary_v_y(&mesh, V2_EXT);
	CustomInitialCondition_e boundary_e(&mesh,KAPPA);

		Solution<double> R_h_1, R_h_2;

 //--------- Filters for visualization of pressure & velocity
 	PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 1);
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 2);

  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView s1("rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  s1.set_min_max_range(0, 1.);
  s2.set_min_max_range(0., 1.);
  s3.set_min_max_range(0., 1.);
	pressure_view.set_min_max_range(0.,1.);	
  if(HERMES_VISUALIZATION){  
			s1.show(&prev_rho);
			s2.show(&vel_x);
			s3.show(&vel_y);
			pressure_view.show(&pressure);
	}

 
	OrderView mview_r("mesh", new WinGeom(0, 0, 500, 400));
	OrderView mview_vx("mesh", new WinGeom(0, 0, 500, 400));
	OrderView mview_vy("mesh", new WinGeom(0, 0, 500, 400));
	OrderView mview_e("mesh", new WinGeom(0, 0, 500, 400));	
  		
	PressureFilter pressure_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), KAPPA);
	VelocityFilter vel_x_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 1);
	VelocityFilter vel_y_low(Hermes::vector<MeshFunction<double>*>(&low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e), 2);
	ScalarView s1_n("low_rho", new WinGeom(0, 0, 600, 300));
	ScalarView s2_n("low_rho_v_x", new WinGeom(700, 0, 600, 300));
	ScalarView s3_n("low_rho_v_y", new WinGeom(0, 400, 600, 300));
	ScalarView s4_n("low_pressure", new WinGeom(700, 400, 600, 300));
	s1_n.set_min_max_range(0., 1.);
	s2_n.set_min_max_range(0., 1.);
	s3_n.set_min_max_range(0., 1.);
	s4_n.set_min_max_range(0., 1.);
//------------


//-----------Weakforms
  EulerEquationsWeakForm_Mass wf_mass(time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K(KAPPA, time_step, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerEquationsWeakForm_K  wf_K_low(KAPPA, time_step, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);
  EulerBoundary wf_boundary(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);
  EulerBoundary wf_boundary_low(KAPPA, &boundary_rho, &boundary_v_x, &boundary_v_y,  &boundary_e, &low_rho, &low_rho_v_x, &low_rho_v_y, &low_rho_e);


  // Set up the solver, matrix, and rhs according to the solver selection.
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ;   
	UMFPackMatrix<double> * matrix_L_low = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * low_matrix = new UMFPackMatrix<double> ;  
	UMFPackMatrix<double> * lowmat_rhs = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double> * matrix_dS_low = new UMFPackMatrix<double> ; 

	double* u_L = NULL; 

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	int ps = 1;
	char title[100]; char filename[40];
	
	AsmList<double>*  al = new AsmList<double>;	
	

//Timestep loop
do
{	 

		info(" Time step %d, time %3.5f, ndof=%i", ts, current_time, ndof); 
		
		
			if ((ts > 1 && ts % UNREF_FREQ == 0)||(space_rho.get_num_dofs() >= NDOF_STOP)) 
    { 
    	info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        case 2: mesh.unrefine_all_elements();
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        case 3: mesh.unrefine_all_elements();
								space_rho.set_mesh(&mesh);			space_rho.set_uniform_order(P_INIT); 
								space_rho_v_x.set_mesh(&mesh);	space_rho_v_x.set_uniform_order(P_INIT);  
								space_rho_v_y.set_mesh(&mesh);	space_rho_v_y.set_uniform_order(P_INIT);  
								space_e.set_mesh(&mesh);	space_e.set_uniform_order(P_INIT);   
                break;
        default: Exceptions::Exception("Wrong global derefinement method.");
      }      
    }	

		ps=1;

	//Adaptivity loop
	do{

		ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
		dof_rho = space_rho.get_num_dofs();
		dof_v_x = space_rho_v_x.get_num_dofs();
		dof_v_y = space_rho_v_y.get_num_dofs();
		dof_e = space_e.get_num_dofs();

  // Initialize the FE problem.
  	DiscreteProblem<double>* dp_mass = new DiscreteProblem<double>(&wf_mass, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
 	 DiscreteProblem<double> *dp_boundary = new DiscreteProblem<double>(&wf_boundary, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
 	 DiscreteProblem<double>* dp_K = new DiscreteProblem<double>(&wf_K, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));


	double* coeff_vec = new double[ndof];
	double* coeff_vec_2 = new double[ndof];
	UMFPackVector<double> * vec_rhs = new UMFPackVector<double> (ndof);
	double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
	double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];	
	double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	


    dp_mass->assemble(mass_matrix);
	//----------------------MassLumping M_L--------------------------------------------------------------------
		UMFPackMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Projection of the initial condition
			Lumped_Projection::project_lumped(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec, matrix_solver);
			
		if(HERMES_VISUALIZATION_LOW)
		{
			Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));	
			s1_n.show(&low_rho);
			s2_n.show(&vel_x_low);
			s3_n.show(&vel_y_low);
			s4_n.show(&pressure_low);
		}					
		OGProjection<double>::project_global(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e),Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);


		dp_boundary->assemble(matrix_dS);
    dp_K->assemble(lowmat_rhs);

					//------------------------artificial DIFFUSION D---------------------------------------		
			UMFPackMatrix<double> * diffusion = artificialDiffusion(KAPPA,coeff_vec,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, lowmat_rhs);
			lowmat_rhs->add_matrix(diffusion); //L(U)
			lowmat_rhs->add_matrix(matrix_dS); //L(U)+dS(U) 

			low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
			low_matrix->multiply_with_Scalar(-theta*time_step);  //-theta L(U)
			low_matrix->add_matrix(lumped_matrix); 				//M_L/t - theta L(U)
			lowmat_rhs->multiply_with_Scalar((1.0-theta)*time_step);  //(1-theta)L(U)
			lowmat_rhs->add_matrix(lumped_matrix);  //M_L/t+(1-theta)L(U)

	//-------------rhs lower Order M_L/tau+ (1-theta)(L) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2); 
			vec_rhs->zero(); vec_rhs->add_vector(coeff_vec_2); 

	//-------------------------solution of lower order------------ (M_L/t - theta L(U))U^L = (M_L/t+(1-theta)L(U))U^n
			UMFPackLinearSolver<double> * lowOrd = new UMFPackLinearSolver<double> (low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_sln_vector();  
        Solution<double>::vector_to_solutions(u_L, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&low_rho,&low_rho_v_x,&low_rho_v_y,&low_rho_e));
			 }else error("Matrix solver failed.\n");


		if(HERMES_VISUALIZATION_LOW){
				s1_n.show(&low_rho);
				s2_n.show(&vel_x_low);
				s3_n.show(&vel_y_low);
				s4_n.show(&pressure_low);
			}


//-------------------------

	if(ps==1){
		HPAdapt * adapting = new HPAdapt(Hermes::vector<Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
		int* elements_to_refine = new int[space_rho.get_mesh()->get_max_element_id()]; 
		int* no_of_refinement_steps = new int[space_rho.get_mesh()->get_max_element_id()];	
		double* elem_error = new double[space_rho.get_mesh()->get_max_element_id()];

		for_all_active_elements(e, space_rho.get_mesh()){	elements_to_refine[e->id] = 2; no_of_refinement_steps[e->id]=0;	}
		calc_elem_error(&space_rho, &low_rho,&R_h_1, &R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error,true);
	//calc_elem_error(&space_rho_v_x, &low_rho_v_x,&R_h_1, &R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error,true);
		//calc_elem_error(&space_rho_v_y, &low_rho_v_y,&R_h_1, &R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error,true);
		//calc_elem_error(&space_e, &low_rho_e,&R_h_1, &R_h_2,adapting,h_min,h_max,elements_to_refine,	no_of_refinement_steps,elem_error,true);

		adapting->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

			delete adapting;
			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
			delete [] elem_error;			

	}else{
	

		//------------as=2-------
		 DiscreteProblem<double>* dp_boundary_low = new DiscreteProblem<double>(&wf_boundary_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
			DiscreteProblem<double>* dp_K_low = new DiscreteProblem<double>(&wf_K_low, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

		
			dp_boundary_low->assemble(matrix_dS_low);	
	 	 	dp_K_low->assemble(matrix_L_low);
			UMFPackMatrix<double> * diffusion_low = artificialDiffusion(KAPPA,u_L,&space_rho,&space_rho_v_x, &space_rho_v_y,&space_e, dof_rho, dof_v_x, dof_v_y,dof_e, matrix_L_low);
			matrix_L_low->add_matrix(diffusion_low); //L(U)
			matrix_L_low->add_matrix(matrix_dS_low); //L(U)+dS(U) 
	
	bool* fct =NULL;	
	
if(FCT_W_HANGING_NODES== false){			
	fct = new bool[ndof];
	for(int i = 0; i< ndof; i++) fct[i] = false;
	p1_list(&space_rho, fct, al,dof_rho, dof_v_x, dof_v_y,dof_e);
	}

			//---------------------------------------antidiffusive fluxes-----------------------------------	
			antidiffusiveFlux(mass_matrix,lumped_matrix,diffusion, matrix_L_low, u_L, coeff_vec_2, P_plus, P_minus, Q_plus, Q_minus,R_plus, R_minus, dof_rho, dof_v_x, dof_v_y,dof_e,fct);
				
			
			for(int i=0; i<ndof;i++){
							 coeff_vec[i] = u_L[i]+ coeff_vec_2[i]*time_step/lumped_matrix->get(i,i);		
			}

			Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&new_rho, &new_rho_v_x, &new_rho_v_y, &new_e));	

 			prev_rho.copy(&new_rho); 
			prev_rho_v_x.copy(&new_rho_v_x); 
			prev_rho_v_y.copy(&new_rho_v_y); 
			prev_e.copy(&new_e); 
			
			
		/*	 prev_rho.copy(&low_rho); 
			prev_rho_v_x.copy(&low_rho_v_x); 
			prev_rho_v_y.copy(&low_rho_v_y); 
			prev_e.copy(&low_rho_e);*/
			
			 // Visualize the solution.
			if(HERMES_VISUALIZATION){
				sprintf(title, "pressure: ts=%i",ts);
				pressure_view.set_title(title);
				s1.show(&prev_rho);
				s2.show(&vel_x);
				s3.show(&vel_y);
				pressure_view.show(&pressure);
  		}
  		
  		      // Output solution in VTK format.
    if((VTK_VISUALIZATION)&&(ts  % VTK_FREQ == 0)) 
    {
    	Linearizer lin_p,lin_rho,lin_v_x; 	Orderizer ord_p;		
      sprintf(filename, "pressure-%i.vtk", ts);
      lin_p.save_solution_vtk(&pressure, filename, "Pressure", true);  
      sprintf(filename, "density-%i.vtk", ts);
      lin_rho.save_solution_vtk(&prev_rho, filename, "density", true); 
      sprintf(filename, "mesh-%i.vtk", ts); 
     	ord_p.save_mesh_vtk(&space_rho, filename);  
         sprintf(filename, "vel_x-%i.vtk", ts); 	
     		lin_v_x.save_solution_vtk(&vel_x,filename, "velocity_x", true);
    }

		delete diffusion_low;
		delete dp_boundary_low;
		delete dp_K_low;
		if(fct!=NULL) delete [] fct;

	
	}




		if(HERMES_VISUALIZATION_ORDER)
		{
			mview_r.show(&space_rho);
			mview_vx.show(&space_rho_v_x);
			mview_vy.show(&space_rho_v_y);
			mview_e.show(&space_e);
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


		ps++;

	}while(ps<3);


	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;




}
while (current_time < T_FINAL);

     if(VTK_VISUALIZATION)
     {
        Linearizer lin_p, lin_v_x, lin_v_y, lin_rho; 	Orderizer ord_p;		
			lin_p.save_solution_vtk(&pressure, "p_end.vtk", "pressure", true);
			lin_v_x.save_solution_vtk(&vel_x, "vx_end.vtk", "velocity_x", true);
			lin_v_y.save_solution_vtk(&vel_y, "vy_end.vtk", "velocity_y",true);
			lin_rho.save_solution_vtk(&prev_rho, "rho_end.vtk", "density", true);
			     ord_p.save_mesh_vtk(&space_rho, "mesh_end.vtk");  
		}


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

