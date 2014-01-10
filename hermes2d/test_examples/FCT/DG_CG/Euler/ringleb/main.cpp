#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "interface.h"
#include "euler_flux.h"

#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"

using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 10;
const double T_FINAL =10000000;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// Kappa.
const double KAPPA = 1.4; 
    
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

bool equal(double x, double y)
{
		if(std::fabs(x-y)< 1e-10) return true;
		else false;
}

int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain2.mesh", basemesh);
Element* e = NULL;
Node* vn=NULL; 



  // Perform initial mesh refinements (optional).
  for (int ref=0; ref < INIT_REF_NUM; ref++)
	{ 
		basemesh->refine_all_elements();
	}


int no_nodes_per_side = std::pow(2,INIT_REF_NUM)+1;
int no_en_per_side = std::pow(2,INIT_REF_NUM);
double h = 2./no_en_per_side;
//printf("nodes=%i", no_nodes_per_side);

Node* node_pair[2][no_en_per_side/2];
double a, rho,q, p, J, k ;
	for_all_vertex_nodes(vn, basemesh)
	{
		if(vn->bnd)
		{
			
			if(std::fabs(vn->y) == 1)
			{
				int int_node= -100;
				for(int i = 0; i< no_nodes_per_side; i++)
				{
						if(equal(vn->x, -1+i*h))
						{	int_node = i;
							break;	
						}
				}
				if(int_node<0) printf("error");
				q= 0.5;
				a = std::sqrt(1.-(KAPPA-1.)/2.*q*q);
				rho = std::pow(a, (2./(KAPPA-1.)));
				p = 1./KAPPA*std::pow(a, (2.*KAPPA/(KAPPA-1.)));
				J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
				k = 1.5-int_node*0.8/no_en_per_side;
				vn->x = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;
				if(vn->y== -1) vn->y= -1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));
				if(vn->y == 1) vn->y= 1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));			

			}else if(vn->y== 0)
			{
				if(equal(vn->x,-1))
				{
					k = 1.5; 
					node_pair[0][0] = vn;					
				}else if(equal(vn->x,1))
				{	
					k = 0.7;
					node_pair[1][0] = vn;				
				}
				q = k;
				a = std::sqrt(1.-(KAPPA-1.)/2.*q*q);
				rho = std::pow(a, (2./(KAPPA-1.)));
				p = 1./KAPPA*std::pow(a, (2.*KAPPA/(KAPPA-1.)));
				J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));			
				vn->x = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;


			}else if(std::fabs(vn->x) == 1.)			
			{
				int int_node= 0; int abs_node = 0;
				for(int i = 1; i< (no_en_per_side)/2.; i++)
				{
						if(equal(vn->y, i*h))
						{	int_node = i;
							abs_node = i;							
							break;	
						}else if(equal(vn->y, -i*h))
						{	int_node = -i;
							abs_node = i;							
							break;	
						}
				}
				if(int_node==0) printf("error");
				if(vn->x<0){k = 1.5; q = k-abs_node*h;  if(int_node>0) node_pair[0][int_node] = vn;
				}else{ k = 0.7; q = k-abs_node*h*0.2;	if(int_node>0) node_pair[1][int_node] = vn; }
				
				if(q<0.5) printf("error");
				a = std::sqrt(1.-(KAPPA-1.)/2.*q*q);
				rho = std::pow(a, (2./(KAPPA-1.)));
				p = 1./KAPPA*std::pow(a, (2.*KAPPA/(KAPPA-1.)));
				J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));			
				vn->x = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;
				if(vn->y<0) vn->y= -1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));
				if(vn->y >0) vn->y= 1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));
			}
		}
	}


Node* vn_1 = NULL; Node* vn_2 = NULL;
int coord_y, coord_x;
	
		for_all_vertex_nodes(vn, basemesh)
		{
			coord_y = 0; coord_x = 0;
			if(!vn->bnd)
			{
				for(int j = 0; j< (no_en_per_side)/2.;j++)
					if(equal(std::fabs(vn->y),j*h)){coord_y = j; break;}
				for(int i = 1; i< (no_en_per_side); i++)
					if(equal(vn->x,-1+i*h)){coord_x = i; break;}


				vn_1 = node_pair[0][coord_y];
				vn_2 = node_pair[1][coord_y];
				double delta_x = (vn_2->x - vn_1->x)/no_en_per_side;
				double delta_y = (vn_2->y - vn_1->y)/no_en_per_side;
				vn->x = vn_1->x + coord_x *delta_x;
				if(vn->y >0) vn->y = vn_1->y + coord_x*delta_y;	
				if(vn->y <0) vn->y = -(vn_1->y +coord_x*delta_y);				

			}

		}

		

				



/* MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(basemesh);
   View::wait();*/

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


  MeshFunctionSharedPtr<double> boundary_rho(new BoundaryCondition_rho(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_x(new  BoundaryCondition_rho_v_x(mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_y(new  ConstantSolution<double>(mesh, V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new BoundaryCondition_rho_e(mesh,KAPPA));	

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

  ScalarView s5("diff_rho", new WinGeom(700, 0, 600, 300));
  ScalarView mach_view("mach", new WinGeom(0, 700, 600, 300));
			//s2.set_min_max_range(-1.,1.);
			//s3.set_min_max_range(-1.,1.);
			//mach_view.set_min_max_range(0.0, 0.1);
			//pressure_view.set_min_max_range(0.68,0.72);
			//s1.set_min_max_range(0.9, 1.6);
s5.set_min_max_range(-0.1,0.1);

                        s1.show(init_rho);
						s4.show(init_e);

Linearizer lin_p, lin_m, lin_rho;



//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
 
	//NumericalFlux* num_flux = new HLLNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new ApproxRoeNumericalFlux(KAPPA, euler_fluxes); 
NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new StegerWarmingNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new VijayasundaramNumericalFlux(KAPPA);
//NumericalFlux* num_flux =new OsherSolomonNumericalFlux(KAPPA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);
	EulerK wf_convection_init(KAPPA,init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

	EulerS wf_bdry_init(KAPPA,mesh, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e,false);
	EulerS wf_bdry(KAPPA,mesh, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,false);





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
double residual = 0.;

		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(mass_matrix);
mass_matrix->multiply_with_Scalar(1./time_step);
CSCMatrix<double> * lumped_matrix = massLumping(mass_matrix);

//Timestep loop
do
{	 
  //if(ts  % 100 == 1)
	//Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, norm =%f, norm_rel = %f", ts, current_time, ndof, norm, norm_rel); 	  
	Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, norm =%f, norm_rel = %f, res = %f", ts, current_time, ndof, norm, norm_rel,residual); 		
 	 if(ts!=1){
			dp.assemble(K_matrix, vec_conv);
			dp_bdry.assemble(dS_matrix, vec_bdry);
			//dp_DG.assemble(vec_dg);	
		  	//dp_DG.assemble(dg_matrix,vec_dg);	
		}else{
			dp_init.assemble(K_matrix, vec_conv);
			dp_bdry_init.assemble(dS_matrix, vec_bdry);
				//dp_DG.assemble(vec_dg);	
		 	//dp_DG_init.assemble(dg_matrix,vec_dg);	
		}

		//(M-theta(K+ds))u(n+1) = Sn +Ku(n) +(M-theta(Kn+ds))u(n)

//matrix->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
//matrix->add_sparse_matrix(K_matrix);

CSCMatrix<double>* diff= NULL;

matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());//L(U) = KU+SU
	diff = artificialDiffusion(KAPPA,coeff_vec,spaces,dof_rho,dof_v_x, dof_v_y, dof_e,K_matrix);
			 matrix->add_sparse_matrix(diff);


		matrix->add_sparse_matrix(dS_matrix);
		matrix->multiply_with_Scalar(-theta);  //-theta L(U)	
		//matrix->add_sparse_matrix(mass_matrix); 			//M/t - theta L(U)
matrix->add_sparse_matrix(lumped_matrix); 


	//-------------rhs: M/tau+ (1-theta)(L) u^n------------		
	//	matrix->multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs->zero(); 
		//vec_rhs->add_vector(coeff_vec_2); 
		//vec_rhs->add_vector(vec_dg); 
		vec_rhs->add_vector(vec_bdry); 
		vec_rhs->add_vector(vec_conv); 
residual = 0;
for(int i = 1; i<ndof; i++)
	residual +=vec_rhs->get(i)*vec_rhs->get(i);


	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
			norm = get_l2_norm(solver->get_sln_vector(), ndof);	
		Solution<double>::vector_to_solutions(solver->get_sln_vector(), spaces, diff_slns);	

		for(int i=0; i<ndof;i++)		
					coeff_vec[i] += solver->get_sln_vector()[i];
							


			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);


	

			// Visualize the solution.
               MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(prev_slns));
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
						s4.show(prev_e);
				s5.show(diff_slns[0]);

	//View::wait(HERMES_WAIT_KEYPRESS);



	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix->free();
if(diff!=NULL) delete diff;
 

double abs = get_l2_norm(coeff_vec, ndof);
norm_rel= norm/abs;


//}while (current_time < T_FINAL);
}while ((current_time < T_FINAL)||(norm <1e-12)||(norm_rel<1e-08));



Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  


		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));


				pressure->reinit();
				mach->reinit();
			  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );       
			lin_p.save_solution_vtk(pressure, filename, "pressure", true);
sprintf(filename, "m-%i.vtk", ts );      
			lin_m.save_solution_vtk(mach, filename, "mach", true);

sprintf(filename, "rho-%i.vtk", ts );  
			lin_rho.save_solution_vtk(prev_slns[0], filename, "density", true);





/*
Orderizer ord_space;
ord_space.save_orders_vtk(spaces[0], "space.vtk");*/



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

