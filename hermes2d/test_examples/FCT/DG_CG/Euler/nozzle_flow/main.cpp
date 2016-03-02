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

const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.
const double time_step = 0.25;
const double T_FINAL = 60000000.;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
const bool DG = true;
     
const bool serendipity= true;
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
  mloader.load("domain2.mesh", basemesh);
Element* e = NULL;Node* vn=NULL;

 //for_all_active_elements(e, basemesh)
	//	basemesh->refine_quad_to_triangles(e);


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();
		//y-Koord. der Knoten auf cos-Kurve verschieben
			 
			for_all_vertex_nodes(vn, basemesh)
			{	
				if(vn->bnd) 
					if((vn->x>0.)&&(vn->x<4.))		
					{		if(vn->y>0.)
							{
								
								vn->y = (Hermes::cos(PI*vn->x/2.)+3.)/4.;
							}else if(vn->y<0.)
									
							vn->y = -(Hermes::cos(PI*vn->x/2.)+3.)/4.;
					}
			}
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
double CFL = (time_step/delta_max*0.2)+1;
printf("CFL = %f \n", CFL);


/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


/*
SpaceSharedPtr<double> space_rho(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new H1Space<double>(mesh, P_INIT));
/*	

	SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));*/
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
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new  CustomInitialCondition_v_x_rho(mesh,KAPPA));	
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
  MeshFunctionSharedPtr<double> boundary_v_x(new  CustomInitialCondition_v_x_rho(mesh,KAPPA));	
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
//ScalarView s4("prev_e", new WinGeom(700, 700, 600, 300));

  ScalarView mach_view("mach", new WinGeom(0, 700, 600, 300));
		/*	s2.set_min_max_range(1.,4.);
			s3.set_min_max_range(-1.,1.);
			mach_view.set_min_max_range(0.17, 0.48);
			pressure_view.set_min_max_range(0.68,0.72);
			s1.set_min_max_range(0.91, 1.);*/


  ScalarView s5("diff_rho", new WinGeom(700, 0, 600, 300));



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
//NumericalFlux* num_flux =new VijayasundaramNumericalFlux(KAPPA);


	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA,mesh, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA,mesh, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);
	EulerK wf_convection_init(KAPPA,init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass;

	EulerS wf_bdry_init(KAPPA,mesh,num_flux, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerS wf_bdry(KAPPA,mesh, num_flux,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);





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
	CSCMatrix<double> * matrix_2 = new CSCMatrix<double>;  

    OGProjection<double> ogProjection;
		SimpleVector<double> * vec_dg = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_res = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_bdry = new SimpleVector<double> (ndof);
		SimpleVector<double> * vec_conv = new SimpleVector<double> (ndof);
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
double norm = 1000;
double norm_rel = 1000;
double residual = 10.;
double residual_norm = 10.;

		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(mass_matrix);
mass_matrix->multiply_with_Scalar(1./time_step);
//CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);
FILE * pFile;
pFile = fopen ("residual.txt","a");
    fprintf (pFile,"DOFS:%i, CFL=%f \n",ndof,CFL);
fclose (pFile);

Linearizer lin;char filename[40];
		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));
		
//---------- for residual-calculation=------------		
ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
 errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(space_rho->get_mesh()));
 Hermes::vector<MeshFunctionSharedPtr<double> > zero_slns(sln_zero, sln_zero, sln_zero, sln_zero);
 //------------------
//Timestep loop
do
{	 
  //if(ts  % 100 == 1)
	//Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, norm =%f, norm_rel = %f, res = %f", ts, current_time, ndof, norm, norm_rel,residual); 
Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, res = %f", ts, current_time, ndof, residual); 	  
 	
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

		if(DG) 
		{
matrix->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
matrix->add_sparse_matrix(K_matrix);
		}else{

matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());//L(U) = KU+SU
		}


	CSCMatrix<double>* diff= NULL;

	/*if(ts<4){
		matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());//L(U) = KU+SU
	diff = artificialDiffusion(KAPPA,coeff_vec,spaces,dof_rho,dof_v_x, dof_v_y, dof_e,K_matrix);
			 matrix->add_sparse_matrix(diff);

		}else{

matrix->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
matrix->add_sparse_matrix(K_matrix);

}*/
		matrix->add_sparse_matrix(dS_matrix);
		matrix->multiply_with_Scalar(-theta);  //-theta L(U)	

/*
if(ts<4)
	{	
//		matrix->add_sparse_matrix(mass_matrix); 			//M/t - theta L(U)
matrix->add_sparse_matrix(lumped_matrix); 

}else */


if((residual_norm<1e-3))
mass_matrix->multiply_with_Scalar(1./10.);

matrix->add_sparse_matrix(mass_matrix); 

/*
//matrix_2->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
matrix_2->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());
matrix_2->zero();
matrix_2->add_sparse_matrix(dS_matrix);
matrix_2->multiply_with_Scalar(-theta);
matrix_2->add_sparse_matrix(mass_matrix); 
*/


	//-------------rhs: M/tau+ (1-theta)(L) u^n------------		
		//matrix_2->multiply_with_vector(coeff_vec, coeff_vec_2);
//diff->multiply_with_vector(coeff_vec, coeff_vec_2);
//mass_matrix->multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs->zero(); 
		//vec_rhs->add_vector(coeff_vec_2); 
		//if(ts>=4) 
	if(DG) vec_rhs->add_vector(vec_dg); 
		vec_rhs->add_vector(vec_bdry); 
		vec_rhs->add_vector(vec_conv); 

vec_res->zero();
if(DG)  vec_res->add_vector(vec_dg); 
		vec_res->add_vector(vec_bdry); 
		vec_res->add_vector(vec_conv); 

//---For residual calculation
for(int i=0; i<ndof;i++)		
					coeff_vec_2[i] =vec_res->get(i);
Solution<double>::vector_to_solutions(coeff_vec_2, spaces, diff_slns);	
errorCalculator_l2.calculate_errors(diff_slns, zero_slns);
double err_l2_2 = errorCalculator_l2.get_total_error_squared();
residual_norm  = Hermes::sqrt(err_l2_2);
//-----------------------------------------------
residual = 0;
for(int i = 1; i<ndof; i++)
	residual +=vec_res->get(i)*vec_res->get(i);

	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
for(int i=0; i<ndof;i++)		
					coeff_vec_2[i]= solver->get_sln_vector()[i];// - coeff_vec[i];
			norm = get_l2_norm(coeff_vec_2, ndof);	
		Solution<double>::vector_to_solutions(coeff_vec_2, spaces, diff_slns);	

		for(int i=0; i<ndof;i++)		
					coeff_vec[i]= solver->get_sln_vector()[i]+ coeff_vec[i];;
							


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

                        MeshFunctionSharedPtr<double> mach(new MachNumberFilter(prev_slns, KAPPA));
                        sprintf(title, "Mach: ts=%i",ts);
                        mach_view.set_title(title);
                        mach->reinit();
                        mach_view.show(mach);

                        sprintf(title, "Density: ts=%i",ts);
                        s1.set_title(title);
                        s1.show(prev_rho);
				s5.show(diff_slns[0]);*/

	//View::wait(HERMES_WAIT_KEYPRESS);

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
/*
if(ts==4)
{ for(int i = 0; i<4; i++) spaces[i]->set_uniform_order(2);
Space<double>::assign_dofs(spaces);
ndof = Space<double>::get_num_dofs(spaces);
			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;
			delete vec_bdry;
			delete vec_conv;

		 vec_dg = new SimpleVector<double> (ndof);
		 vec_rhs = new SimpleVector<double> (ndof);
		 vec_bdry = new SimpleVector<double> (ndof);
		 vec_conv = new SimpleVector<double> (ndof);
		 coeff_vec = new double[ndof];	
		 coeff_vec_2 = new double[ndof];

  ogProjection.project_global(spaces,prev_slns, coeff_vec, norms_l2 );
		Space<double>::assign_dofs(spaces);
			dp_mass.set_spaces(spaces);
			dp.set_spaces(spaces);
			dp_bdry.set_spaces(spaces);
			dp_DG.set_spaces(spaces);
		  dp_mass.assemble(mass_matrix);

//mass_matrix->multiply_with_Scalar(1./1e-2);

}*/

		delete solver;
		matrix->free();
matrix_2->free();
if(diff!=NULL) delete diff;
 

double abs = get_l2_norm(coeff_vec, ndof);
norm_rel= norm/abs;

int bound = 0;
for(int i = 0;	i<10; i++)
{
	if(residual < 1./std::pow(10,i)) bound = i;
	else break;
}
Hermes::Mixins::Loggable::Static::info("res = %f, %f < 10^(-%i)", residual, residual_norm, bound); 	  
 	

pFile = fopen ("residual_norm.txt","a");
   // fprintf (pFile,"%i: res = %e < 10^(-%i),residual_norm=%e,  norm =%e, norm_rel = %e \n",ts, residual, bound,residual_norm, norm, norm_rel);
	     fprintf (pFile,"%i: %e\n",ts,residual_norm);
fclose (pFile);

pFile = fopen ("residual.txt","a");
    fprintf (pFile,"%i: res = %e < 10^(-%i),residual_norm=%e, \n",ts, residual, bound,residual_norm);
fclose (pFile);

if(ts==100){
			pressure->reinit();
				mach->reinit();
			  sprintf(filename, "p-%i.vtk", ts );        
			lin.save_solution_vtk(pressure, filename, "pressure", false);
sprintf(filename, "m-%i.vtk", ts );
			lin.save_solution_vtk(mach, filename, "mach", false);
sprintf(filename, "rho-%i.vtk", ts );       
			lin.save_solution_vtk(prev_slns[0], filename, "density", false);

}

}
//while (ts < 20);
//while ((current_time < T_FINAL)&&(residual>1e-12));
//while ((ts < 50)&&(residual>1e-10));
while ((ts <200)&&(residual_norm>1e-10));

//if(residual<=1e-8) printf("Residual small enough");


//Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

	//	MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		//MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));

				pressure->reinit();
				mach->reinit();

        Linearizer lin_p;
			lin_p.save_solution_vtk(pressure, "p_end.vtk", "pressure", false);
        Linearizer lin_m;
			lin_m.save_solution_vtk(mach, "m_end.vtk", "mach",false);
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_slns[0], "rho_end.vtk", "density", false);



	/*		  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );

        Linearizer lin_p;
			lin_p.save_solution_vtk(pressure, filename, "pressure", true);
sprintf(filename, "m-%i.vtk", ts );

        Linearizer lin_m;
			lin_m.save_solution_vtk(mach, filename, "mach", true);

sprintf(filename, "rho-%i.vtk", ts );
        Linearizer lin_rho;
			lin_rho.save_solution_vtk(prev_slns[0], filename, "density", true);*/


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

