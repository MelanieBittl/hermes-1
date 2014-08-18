#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "interface.h"
#include "euler_flux.h"
#include "lumped_projection.h"



using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"

const int INIT_REF_NUM =1;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 3.65;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// GAMMA.
const double GAMMA = 1.4; 

const double R = 287;

// Penalty Parameter.
double SIGMA = 10000000.;

//Particle density_particle
const double density_particle = 10.; 

bool view_3D = false;

const double particle_start = 0.2;

    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

    static int verf_crit(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
        if(e->vn[i]->x >29.14)
          return 2; //vertikale verfeinerung
      return -1;
    }	
        static int verf_crit_2(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
		{ if((e->vn[i]->y >=15.7)||(e->vn[i]->y <12))
          return 0;
		  else return 2;
		}
      return -1;
    }	
    
            static int verf_crit_3(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
		{ if((e->vn[i]->y >=15.7)&&(e->vn[i]->y<31.35)||(e->vn[i]->y <12)&&(e->vn[i]->y>=-3.65))
	          return 0;
		  else return 2;
		}
      return -1;
    }	
    
    static int verf_crit_4(Element* e)
    {
		 bool all = false;
      for (unsigned int i = 0; i < e->get_nvert(); i++)
		{ if((e->vn[i]->y <=15.35)&&(e->vn[i]->y >=12.35))
          all = true;
			else return -1;
		}
		if(all ==true) return 1;
      return -1;
    }	
    
    
   static int verf_crit_5(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
		{ if((e->vn[i]->y >=15.35)&&(e->vn[i]->y<31.35)||(e->vn[i]->y <12)&&(e->vn[i]->y>=-3.65))
          return 0;
		  else return 2;
		}
      return -1;
    }	
    
       static int verf_crit_6(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
		{ if((e->vn[i]->x >21.028)||(e->vn[i]->x<13.62))
          return 0;
		  else return 2;
		}
      return -1;
    }	
    
int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
 //mloader.load("domain_all.mesh", basemesh);
  mloader.load("mesh_new60.mesh", basemesh);

  // Perform initial mesh refinements (optional).
Element* e = NULL;Node* vn=NULL;
//basemesh->refine_by_criterion(verf_crit,3);

  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
		//basemesh->refine_all_elements();

		}
	
basemesh->refine_by_criterion(verf_crit,3);
basemesh->refine_by_criterion(verf_crit_5,2);

//basemesh->refine_all_elements();



//basemesh->refine_by_criterion(verf_crit_2,2);
//basemesh->refine_by_criterion(verf_crit_3,1);
//basemesh->refine_by_criterion(verf_crit_4,1);



 	 mesh->copy(basemesh);

Element* test_element = RefMap::element_on_physical_coordinates(true, mesh, 16, 14.);






		SpaceSharedPtr<double> space_rho_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_g(new H1Space<double>(mesh, P_INIT));


		SpaceSharedPtr<double> space_rho_p(new H1Space<double>(mesh, P_INIT));	
//Orderizer ord;
//ord.save_mesh_vtk(space_rho_p, "mesh.vtk");
/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();
*/

int dof_rho = space_rho_g->get_num_dofs();
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho_g, space_rho_v_x_g, space_rho_v_y_g, space_e_g,space_rho_p);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho_g(new CustomInitialCondition_rho(mesh,  GAMMA));	
  MeshFunctionSharedPtr<double> init_rho_v_x_g(new  CustomInitialCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_v_y_g(new  CustomInitialCondition_rho_v_y(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_e_g(new CustomInitialCondition_e(mesh, GAMMA));

 MeshFunctionSharedPtr<double> init_rho_p(new ConstantSolution<double>(mesh, 0.));	


  MeshFunctionSharedPtr<double> bdry_rho_g(new BoundaryCondition_rho(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new  BoundaryCondition_rho_v_y(mesh, GAMMA));	
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
  ScalarView temp_view_g("Temperature-gas", new WinGeom(200, 750, 400, 300));

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

  EulerSource wf_source_init(density_particle, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p);
EulerSource wf_source(density_particle,prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p);
EulerSource wf_source_low(density_particle,low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g,low_rho_p);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);

  DiscreteProblem<double> dp_conv_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp_conv(&wf_convection, spaces); 

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 


  DiscreteProblem<double> dp_source_init(&wf_source_init,spaces); 
  DiscreteProblem<double> dp_source(&wf_source, spaces); 

DiscreteProblem<double> dp_source_low(&wf_source_low, spaces); 
DiscreteProblem<double> dp_bdry_low(&wf_bdry_low, spaces); 



  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;
		CSCMatrix<double> matrix2;    
	CSCMatrix<double>  mass_matrix;   
	CSCMatrix<double>  source_matrix;  
	CSCMatrix<double>  conv_matrix; 
	CSCMatrix<double>  bdry_matrix;  
CSCMatrix<double> matrixL; 

    OGProjection<double> ogProjection;
Lumped_Projection lumpedProjection;
	
		SimpleVector<double>  vec_rhs(ndof);
		SimpleVector<double>  vec_bdry(ndof);
		SimpleVector<double>  vec_conv(ndof);
		SimpleVector<double>  vec_source(ndof);
SimpleVector<double>  vec_res(ndof);
		vec_source.zero(); 

		double* coeff_vec= new double[ndof];	
		double* coeff_vec_2 = new double[ndof];
		double* u_L= new double[ndof];;
		double* flux_vec = new double[ndof];
			double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
			double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];
			double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	



//Projection of the initial condition
 // ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
lumpedProjection.project_lumped(spaces, init_slns, coeff_vec);
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);






		MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach_g(new  MachNumberFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
		MeshFunctionSharedPtr<double> vel_y(new VelocityFilter_y(prev_slns));
		MeshFunctionSharedPtr<double> temp_g(new TempFilter(prev_slns,R, GAMMA));
Linearizer lin;	

	
char title[100];

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
		
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
 	// if(ts!=1){
			dp_conv.assemble(&conv_matrix, &vec_conv);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);
			//dp_source.assemble(&source_matrix, &vec_source);
				if(mach>=particle_start) dp_source.assemble(&vec_source);
		
		/*}else{
			dp_conv_init.assemble(&conv_matrix, &vec_conv);
			dp_bdry_init.assemble(&bdry_matrix, &vec_bdry);
				if(mach>=0.9) dp_source.assemble(&vec_source);
			//dp_source_init.assemble(&source_matrix, &vec_source);
	
		}*/

CSCMatrix<double>* diff = NULL;	
 diff =  artificialDiffusion(GAMMA,coeff_vec,spaces,&conv_matrix);

matrixL.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());
matrixL.add_sparse_matrix(diff);



matrix.create(matrixL.get_size(),matrixL.get_nnz(), matrixL.get_Ap(), matrixL.get_Ai(),matrixL.get_Ax());
		matrix.add_sparse_matrix(&bdry_matrix);	
		matrix.multiply_with_Scalar(-theta); 
		matrix.add_sparse_matrix(lumped_matrix);
 
		//matrix.add_sparse_matrix(&mass_matrix); 


matrix2.create(matrix.get_size(),matrix.get_nnz(), matrix.get_Ap(), matrix.get_Ai(),matrix.get_Ax());	
matrix2.add_sparse_matrix(diff);

	//-------------rhs: ------------		
		matrix2.multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs.zero(); 
		vec_rhs.add_vector(coeff_vec_2); 		
		vec_rhs.add_vector(&vec_bdry); 
		vec_rhs.add_vector(&vec_conv); 
		if(mach>=particle_start) 
		{vec_rhs.add_vector(&vec_source); 
			
		//for(int i = 1; i<ndof; i++)
	//if(vec_source.get(i)<0) printf("source = %f", vec_source.get(i));
		}


vec_res.zero(); vec_res.add_vector(&vec_bdry); vec_res.add_vector(&vec_conv); 
residual = 0;
for(int i = 1; i<ndof; i++)
	residual +=vec_res.get(i)*vec_res.get(i);

	//-------------------------solution of (M-theta(K+B+S)) (u(n+1)-u(n) = Sn +Ku(n) +Bn------------ 
		// Hermes::Mixins::Loggable::Static::info("Solving"); 	
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (&matrix,&vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	

for(int i=0; i<ndof;i++) //coeff_vec[i]= solver->get_sln_vector()[i];
	u_L[i] = solver->get_sln_vector()[i];



Solution<double>::vector_to_solutions(u_L, spaces, low_slns);
			dp_bdry_low.assemble(&vec_bdry);			

antidiffusiveFlux(&mass_matrix,lumped_matrix,diff,&matrixL, &vec_bdry,&vec_source, u_L, flux_vec, P_plus,P_minus, Q_plus, Q_minus,  R_plus, R_minus,dof_rho,time_step,GAMMA );

		for(int i=0; i<ndof;i++)	
			coeff_vec[i] = u_L[i]+ flux_vec[i]/lumped_matrix->get(i,i);
				
					

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

if(mach<particle_start)
{	double rho		= prev_rho_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0);  
	double rho_v_x	= prev_rho_v_x_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	double rho_v_y	= prev_rho_v_y_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	double rho_e	= prev_rho_e_g.get_solution()->get_ref_value(test_element, 1,0, 0, 0); 
	mach = QuantityCalculator::calc_mach(rho, rho_v_x, rho_v_y, rho_e, GAMMA);
	//printf("mach = %f", mach);
}
//if((ts>100)&&(mach < particle_start)) mach = 1.;
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
                        s1_g.show(prev_rho_g);
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
								 
								 						sprintf(title, "density_particle: ts=%i",ts);
                        s1_p.set_title(title);
                        s1_p.show(prev_rho_p);
								
*/


	//View::wait(HERMES_WAIT_KEYPRESS);




	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix.free();
		matrix2.free();
		matrixL.free();
		if(diff!=NULL) delete diff;
	if(ts%2000 ==1)
{char filename[40];
		/*pressure_g->reinit();
		mach_g->reinit();
		vel_x->reinit();
		vel_y->reinit();
		temp_g->reinit();		
		sprintf(filename, "p-%i.vtk", ts );       
		lin.save_solution_vtk(pressure_g, filename, "pressure", view_3D);
		sprintf(filename, "m-%i.vtk", ts );      
		lin.save_solution_vtk(mach_g, filename, "mach", view_3D);
		sprintf(filename, "vel_x-%i.vtk", ts );  
		lin.save_solution_vtk(vel_x, filename, "vel_x", view_3D);
		sprintf(filename, "vel_y-%i.vtk", ts );  
		lin.save_solution_vtk(vel_y, filename, "vel_y", view_3D);
		sprintf(filename, "rho_g-%i.vtk", ts );  
		lin.save_solution_vtk(prev_slns[0], filename, "density_gas", view_3D);*/
		sprintf(filename, "rho_p-%i.vtk", ts );  
		lin.save_solution_vtk(prev_slns[4], filename, "density_particle", view_3D);
		//sprintf(filename, "temp-%i.vtk", ts - 1);
		//lin.save_solution_vtk(temp_g, filename, "temp", view_3D);

}

 



}while ((current_time < T_FINAL)&&(residual > 1e-8));




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

				pressure_g->reinit();
				mach_g->reinit();
				vel_x->reinit();
				vel_y->reinit();
				temp_g->reinit();
 
			lin.save_solution_vtk(pressure_g, "pressure_end.vtk", "pressure", view_3D);  
			lin.save_solution_vtk(mach_g, "mach_end.vtk", "mach", view_3D); 
			lin.save_solution_vtk(prev_slns[0], "rho_g_end.vtk", "density_gas", view_3D);  
			lin.save_solution_vtk(prev_slns[4], "rho_p_end.vtk", "density_particle", view_3D);
			lin.save_solution_vtk(vel_x, "vx_end.vtk", "velocity_x", view_3D);
			lin.save_solution_vtk(vel_y, "vy_end.vtk", "velocity_y",view_3D);
			lin.save_solution_vtk(temp_g, "temp_end.vtk", "temp",view_3D);




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

