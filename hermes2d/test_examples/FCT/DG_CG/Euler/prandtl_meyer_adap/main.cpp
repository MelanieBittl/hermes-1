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
const int P_INIT = 1;       						// Initial polynomial degree.
const double time_step = 0.001;
const double T_FINAL = 30000000.;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 const bool DG = true;
     
const bool serendipity =true;


// Adaptivity
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_METHOD = 1;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.2;                      // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const CandList CAND_LIST = H2D_H_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 3;												// max. numbers of adaptivity steps




// Kappa.
const double KAPPA = 1.4; 
// Inlet x-velocity (dimensionless).
const double V1_EXT =3.5;        
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
 // mloader.load("quad.mesh", basemesh);
mloader.load("domain.mesh", basemesh);
  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements();

Element* e;

 	 mesh->copy(basemesh);
	 
	 MeshSharedPtr rho_mesh(new Mesh);//, v_x_mesh(new Mesh), v_y_mesh(new Mesh), e_mesh(new Mesh);
	 
	rho_mesh->copy(basemesh);//v_x_mesh->copy(basemesh); v_y_mesh->copy(basemesh); e_mesh->copy(basemesh); 

double delta_x = 100; double delta_max= 0.; 
for_all_active_elements(e, basemesh)
{
		for(int i = 1; i<e->get_nvert(); i++)
		{
			delta_x = std::fabs(e->vn[i]->x - e->vn[i-1]->x);
			if(delta_x>delta_max) delta_max = delta_x;
		}
}
double CFL = (time_step/delta_max*2.5)+1;
printf("CFL = %f \n", CFL);




SpaceSharedPtr<double> space_rho(new L2_SEMI_CG_Space<double>(rho_mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(rho_mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(rho_mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(rho_mesh, P_INIT, serendipity));
/*SpaceSharedPtr<double> space_rho_v_x(new L2_SEMI_CG_Space<double>(v_x_mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y(new L2_SEMI_CG_Space<double>(v_y_mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e(new L2_SEMI_CG_Space<double>(e_mesh, P_INIT, serendipity));*/

		SpaceSharedPtr<double> space_rho_h1(new H1Space<double>(basemesh, 1));	
		SpaceSharedPtr<double> space_rho_v_x_h1(new H1Space<double>(basemesh,1));	
		SpaceSharedPtr<double> space_rho_v_y_h1(new H1Space<double>(basemesh, 1));	
		SpaceSharedPtr<double> space_e_h1(new H1Space<double>(basemesh, 1));

	/*int dof_rho = space_rho->get_num_dofs();
	int dof_v_x = space_rho_v_x->get_num_dofs();
	int dof_v_y = space_rho_v_y->get_num_dofs();
	int dof_e = space_e->get_num_dofs();*/

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
/* MeshFunctionSharedPtr<double> init_rho(new CustomInitialCondition_rho(rho_mesh,KAPPA));	
  MeshFunctionSharedPtr<double> init_rho_v_x(new  ConstantSolution<double>(v_x_mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> init_rho_v_y(new  ConstantSolution<double>(v_y_mesh,  V2_EXT));	
  MeshFunctionSharedPtr<double> init_e(new CustomInitialCondition_e(e_mesh,KAPPA));*/

  	MeshFunctionSharedPtr<double> init_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> init_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> init_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> init_e(new Solution<double>); 
  dynamic_cast<Solution<double>* >(init_rho.get())->load("rho.xml",space_rho_h1); 
dynamic_cast<Solution<double>* >(init_rho_v_x.get())->load("v_x.xml",space_rho_v_x_h1); 
dynamic_cast<Solution<double>* >(init_rho_v_y.get())->load("v_y.xml",space_rho_v_y_h1); 
dynamic_cast<Solution<double>* >(init_e.get())->load("energy.xml",space_e_h1);


  	MeshFunctionSharedPtr<double> prev_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_e(new Solution<double>);
	 
	  MeshFunctionSharedPtr<double> ref_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> ref_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> ref_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> ref_e(new Solution<double>);
	 
	 	  MeshFunctionSharedPtr<double> sln_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> sln_e(new Solution<double>);


  MeshFunctionSharedPtr<double> boundary_rho(new CustomInitialCondition_rho(rho_mesh,KAPPA));	
  MeshFunctionSharedPtr<double> boundary_v_x(new  ConstantSolution<double>(rho_mesh,  V1_EXT));	
  MeshFunctionSharedPtr<double> boundary_v_y(new  ConstantSolution<double>(rho_mesh, V2_EXT));	
  MeshFunctionSharedPtr<double> boundary_e(new CustomInitialCondition_e(rho_mesh,KAPPA));	

  	MeshFunctionSharedPtr<double> diff_rho(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_rho_v_x(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_rho_v_y(new Solution<double>);
    MeshFunctionSharedPtr<double> diff_e(new Solution<double>);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
		Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);
		Hermes::vector<MeshFunctionSharedPtr<double> > ref_slns(ref_rho, ref_rho_v_x, ref_rho_v_y, ref_e);
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
mach_view.set_min_max_range(2.5, 3.26);
		/*	s2.set_min_max_range(1.,4.);
			s3.set_min_max_range(-1.,1.);
			
			pressure_view.set_min_max_range(0.68,0.72);
			s1.set_min_max_range(0.91, 1.);*/

		OrderView ref_mview("ref_mesh", new WinGeom(500, 0, 500, 400));
  ScalarView s5("diff_rho", new WinGeom(700, 0, 600, 300));

/*s1.show(init_rho);
s2.show(init_rho_v_x);
s3.show(init_rho_v_y);*/

//------------

//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(KAPPA);
NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(KAPPA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(KAPPA);

	EulerInterface wf_DG_init(KAPPA,mesh, init_rho, init_rho_v_x, init_rho_v_y, init_e,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(KAPPA,mesh, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e,num_flux,euler_fluxes,riemann_invariants);
	EulerK wf_convection_init(KAPPA,init_rho, init_rho_v_x, init_rho_v_y, init_e);
	EulerK wf_convection(KAPPA, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  EulerEquationsWeakForm_Mass wf_mass;

	EulerS wf_bdry_init(KAPPA,mesh,num_flux, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e, false);
	EulerS wf_bdry(KAPPA,mesh, num_flux,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, false);
  
  //	EulerS wf_bdry_init(KAPPA,mesh,num_flux, boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, init_rho, init_rho_v_x, init_rho_v_y, init_e, true);
	//EulerS wf_bdry(KAPPA,mesh, num_flux,boundary_rho, boundary_v_x, boundary_v_y,  boundary_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, true);



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

		Space<double>::assign_dofs(spaces);


FILE * pFile;
pFile = fopen ("residual.txt","a");
    fprintf (pFile,"DOFS:%i, CFL=%f \n",ndof,CFL);
fclose (pFile);
//---------- for residual-calculation=------------		
double residual_norm = 10.;
ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
 errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));

 //------------------

	  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);
       //selector.set_error_weights(1.0,1.0,1.0); 	
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);  
  Adapt<double> adaptivity(spaces, &errorCalculator);
  AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
  adaptivity.set_strategy(&stoppingCriterion);
//Timestep loop
do
{	 

Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i, res = %e", ts, current_time, ndof, residual); 
    // Periodic global derefinement. 
   if ((ts > 1 && ts % UNREF_FREQ == 0)||(Space<double>::get_num_dofs(spaces) >= NDOF_STOP)) 
    { 	rho_mesh->copy(basemesh);//v_x_mesh->copy(basemesh); v_y_mesh->copy(basemesh); e_mesh->copy(basemesh); 
		 for(int i=0;i<4;i++) spaces[i]->set_uniform_order(P_INIT);     
    Space<double>::assign_dofs(spaces);	      
    }
    
    
        bool done = false; int as = 1;
	do 
    	{
			
	 Mesh::ReferenceMeshCreator rho_ref_mesh_creator(rho_mesh);
    MeshSharedPtr rho_ref_mesh = rho_ref_mesh_creator.create_ref_mesh();	
    Space<double>::ReferenceSpaceCreator rho_ref_space_creator(space_rho, rho_ref_mesh);
    SpaceSharedPtr<double> rho_ref_space = rho_ref_space_creator.create_ref_space();
	
		    Space<double>::ReferenceSpaceCreator rho_v_x_ref_space_creator(space_rho_v_x, rho_ref_mesh);
    SpaceSharedPtr<double> rho_v_x_ref_space = rho_v_x_ref_space_creator.create_ref_space();
			    Space<double>::ReferenceSpaceCreator rho_v_y_ref_space_creator(space_rho_v_y, rho_ref_mesh);
    SpaceSharedPtr<double> rho_v_y_ref_space = rho_v_y_ref_space_creator.create_ref_space();
	 			    Space<double>::ReferenceSpaceCreator rho_e_ref_space_creator(space_e, rho_ref_mesh);
    SpaceSharedPtr<double> rho_e_ref_space = rho_e_ref_space_creator.create_ref_space();
	
		/* Mesh::ReferenceMeshCreator x_ref_mesh_creator(v_x_mesh);
    MeshSharedPtr x_ref_mesh = x_ref_mesh_creator.create_ref_mesh();
	    Space<double>::ReferenceSpaceCreator rho_v_x_ref_space_creator(space_rho, x_ref_mesh);
    SpaceSharedPtr<double> rho_v_x_ref_space = rho_v_x_ref_space_creator.create_ref_space();
	
			 Mesh::ReferenceMeshCreator y_ref_mesh_creator(v_y_mesh);
    MeshSharedPtr y_ref_mesh = y_ref_mesh_creator.create_ref_mesh();
		    Space<double>::ReferenceSpaceCreator rho_v_y_ref_space_creator(space_rho, y_ref_mesh);
    SpaceSharedPtr<double> rho_v_y_ref_space = rho_v_y_ref_space_creator.create_ref_space();
	 
			 Mesh::ReferenceMeshCreator e_ref_mesh_creator(e_mesh);
    MeshSharedPtr e_ref_mesh = e_ref_mesh_creator.create_ref_mesh(); 
	Space<double>::ReferenceSpaceCreator rho_e_ref_space_creator(space_rho, e_ref_mesh);
    SpaceSharedPtr<double> rho_e_ref_space = rho_e_ref_space_creator.create_ref_space();*/

	/*int ref_dof_rho = rho_ref_space->get_num_dofs();
	int ref_dof_dof_v_x = rho_v_x_ref_space->get_num_dofs();
	int ref_dof_dof_v_y = rho_v_y_ref_space->get_num_dofs();
	int ref_dof_dof_e = rho_e_ref_space->get_num_dofs();	*/	

	  Hermes::vector<SpaceSharedPtr<double> > ref_spaces(rho_ref_space, rho_v_x_ref_space, rho_v_y_ref_space, rho_e_ref_space);
	Space<double>::assign_dofs(ref_spaces);
   int ref_ndof = Space<double>::get_num_dofs(ref_spaces);
  Hermes::Mixins::Loggable::Static::info("ref_ndof: %d \n", ref_ndof);
  ref_mview.show(ref_spaces[3]);
  
  		SimpleVector<double> * vec_dg = new SimpleVector<double> (ref_ndof);
		SimpleVector<double> * vec_rhs = new SimpleVector<double> (ref_ndof);
		SimpleVector<double> * vec_res = new SimpleVector<double> (ref_ndof);
		SimpleVector<double> * vec_bdry = new SimpleVector<double> (ref_ndof);
		SimpleVector<double> * vec_conv = new SimpleVector<double> (ref_ndof);
		double* coeff_vec = new double[ref_ndof];	
		double* coeff_vec_2 = new double[ref_ndof];
		
	if(ts==1){
		  ogProjection.project_global(ref_spaces,init_slns, coeff_vec, norms_l2 );
		  	Space<double>::assign_dofs(ref_spaces);
			Solution<double>::vector_to_solutions(coeff_vec, ref_spaces, prev_slns);			
			Space<double>::assign_dofs(ref_spaces);
		
	}
 	
 	// if(ts!=1){
dp.set_spaces(ref_spaces); 
dp_bdry.set_spaces(ref_spaces); 
			dp.assemble(K_matrix, vec_conv);
			dp_bdry.assemble(dS_matrix, vec_bdry);
			if(DG)
			{dp_DG.set_spaces(ref_spaces); 
				dp_DG.assemble(dg_matrix,vec_dg);	
			}
	/*	}else{
			//Projection of the initial condition
  ogProjection.project_global(ref_spaces,init_slns, coeff_vec, norms_l2 );
			Solution<double>::vector_to_solutions(coeff_vec, ref_spaces, prev_slns);
			
			dp_init.set_spaces(ref_spaces); 
			dp_bdry_init.set_spaces(ref_spaces); 
			dp_init.assemble(K_matrix, vec_conv);
			dp_bdry_init.assemble(dS_matrix, vec_bdry);
		 	if(DG)
			{	dp_DG_init.set_spaces(ref_spaces); 
				dp_DG_init.assemble(dg_matrix,vec_dg);
			}
		}*/
		
		dp_mass.set_spaces(ref_spaces);
		  dp_mass.assemble(mass_matrix);
		mass_matrix->multiply_with_Scalar(1./time_step);
		
		
		//(M-theta(K+ds))u(n+1) = Sn +Ku(n) +(M-theta(Kn+ds))u(n)
  //Hermes::Mixins::Loggable::Static::info("step1 \n");
		if(DG) {
			matrix->create(dg_matrix->get_size(),dg_matrix->get_nnz(), dg_matrix->get_Ap(), dg_matrix->get_Ai(),dg_matrix->get_Ax());
			matrix->add_sparse_matrix(K_matrix);
		}else{
			matrix->create(K_matrix->get_size(),K_matrix->get_nnz(), K_matrix->get_Ap(), K_matrix->get_Ai(),K_matrix->get_Ax());//L(U) = KU+SU
		}

  //Hermes::Mixins::Loggable::Static::info("step2 \n");
		matrix->add_sparse_matrix(dS_matrix);
		matrix->multiply_with_Scalar(-theta);  //-theta L(U)	

//if((residual_norm<1e-3))
//mass_matrix->multiply_with_Scalar(1./10.);

matrix->add_sparse_matrix(mass_matrix); 

  //Hermes::Mixins::Loggable::Static::info("step3 \n");
	//-------------rhs: M/tau+ (1-theta)(L) u^n------------		
		vec_rhs->zero(); 
		if(DG)  vec_rhs->add_vector(vec_dg); 
		vec_rhs->add_vector(vec_bdry); 
		vec_rhs->add_vector(vec_conv); 

vec_res->zero();
vec_res->add_vector(vec_rhs); 

//---For residual calculation
for(int i=0; i<ref_ndof;i++)		
					coeff_vec_2[i] =vec_res->get(i);
Solution<double>::vector_to_solutions(coeff_vec_2, ref_spaces, diff_slns);	
 MeshFunctionSharedPtr<double> sln_zero1(new ZeroSolution<double>(rho_ref_mesh));
  /*MeshFunctionSharedPtr<double> sln_zero2(new ZeroSolution<double>(x_ref_mesh));
   MeshFunctionSharedPtr<double> sln_zero3(new ZeroSolution<double>(y_ref_mesh));
	 MeshFunctionSharedPtr<double> sln_zero4(new ZeroSolution<double>(e_ref_mesh));*/
 Hermes::vector<MeshFunctionSharedPtr<double> > zero_slns(sln_zero1, sln_zero1, sln_zero1, sln_zero1);

errorCalculator_l2.calculate_errors(diff_slns, zero_slns);
double err_l2_2 = errorCalculator_l2.get_total_error_squared();
residual_norm  = Hermes::sqrt(err_l2_2);
//-----------------------------------------------
residual = 0;
for(int i = 0; i<ref_ndof; i++)
	residual +=vec_res->get(i)*vec_res->get(i);

	//-------------------------solution of lower order------------ (M/t - theta L(U))U^L = (M/t+(1-theta)L(U))U^n
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (matrix,vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
for(int i=0; i<ref_ndof;i++)		
					coeff_vec_2[i]= solver->get_sln_vector()[i]; //- coeff_vec[i];
			norm = get_l2_norm(coeff_vec_2, ref_ndof);	
		Solution<double>::vector_to_solutions(coeff_vec_2, ref_spaces, diff_slns);	

		for(int i=0; i<ref_ndof;i++)		
					coeff_vec[i]+= solver->get_sln_vector()[i];						


			Solution<double>::vector_to_solutions(coeff_vec, ref_spaces, ref_slns);
		
		     ogProjection.project_global(ref_spaces, ref_slns, slns);
			  
			      errorCalculator.calculate_errors(slns, ref_slns, true);
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    // If err_est too large, adapt the mesh->
  if((err_est_rel_total < ERR_STOP)||(as>=ADAPSTEP_MAX)){
      done = true;as=1;
	 Solution<double>::vector_to_solutions(coeff_vec, ref_spaces, prev_slns);
  }else
    {as++;
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      Hermes::vector<RefinementSelectors::Selector<double> *> selectors(&selector, &selector,&selector,&selector);
      done = adaptivity.adapt(selectors);
    }

	

			// Visualize the solution.
               /*MeshFunctionSharedPtr<double> vel_x(new VelocityFilter_x(ref_slns));
                MeshFunctionSharedPtr<double> vel_y (new VelocityFilter_y(ref_slns));
                        sprintf(title, "vx: ts=%i",ts);        
                        s2.set_title(title);
                        sprintf(title, "vy: ts=%i",ts);        
                        s3.set_title(title);
                        vel_x->reinit();
                        vel_y->reinit();
                        s2.show(vel_x);
                        s3.show(vel_y);
                       MeshFunctionSharedPtr<double> pressure(new PressureFilter(ref_slns, KAPPA));
                        sprintf(title, "Pressure: ts=%i",ts);
                        pressure_view.set_title(title);
                        pressure->reinit();
                        pressure_view.show(pressure);*/

                        MeshFunctionSharedPtr<double> mach(new MachNumberFilter(ref_slns, KAPPA));
                        sprintf(title, "Mach: ts=%i",ts);
                        mach_view.set_title(title);
                        mach->reinit();
                        mach_view.show(mach);

                        sprintf(title, "Density: ts=%i",ts);
                        s1.set_title(title);
                        s1.show(ref_rho);
				//s5.show(diff_slns[0]);

	//View::wait(HERMES_WAIT_KEYPRESS);
				
			delete solver;
			matrix->free();
			matrix_2->free();
		

			delete[] coeff_vec_2;
			delete [] coeff_vec;
			delete vec_rhs;
			delete vec_dg;
			delete vec_bdry;
			delete vec_conv;
			delete vec_res;

			
		}while (done == false);

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;



 



int bound = 0;
for(int i = 0;	i<10; i++)
{
	if(residual_norm < 1./std::pow(10,i)) bound = i;
	else break;
}
Hermes::Mixins::Loggable::Static::info("res_norm = %e < 10^(-%i)", residual_norm, bound); 	  
 	

pFile = fopen ("residual_norm.txt","a");
	     fprintf (pFile,"%i: %e\n",ts,residual_norm);
fclose (pFile);

pFile = fopen ("residual.txt","a");
    fprintf (pFile,"%i: res = %e < 10^(-%i),residual_norm=%e, \n",ts, residual, bound,residual_norm);
fclose (pFile);



}//while ((current_time < T_FINAL)||(  norm <1e-12)||(norm_rel<1e-08));
//while ((ts < 1000)&&(residual_norm>1e-10));
while (ts < 100);


if(residual<=1e-8) printf("Residual small enough");


Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, KAPPA));
Linearizer lin;

				pressure->reinit();
				mach->reinit();

       
			lin.save_solution_vtk(pressure, "p_end2d.vtk", "pressure", false);    
			lin.save_solution_vtk(mach, "m_end2d.vtk", "mach", false);
			lin.save_solution_vtk(prev_slns[0], "rho_end2d.vtk", "density", false);
			
		lin.save_solution_vtk(pressure, "p_end3d.vtk", "pressure", true);    
			lin.save_solution_vtk(mach, "m_end3d.vtk", "mach", true);
			lin.save_solution_vtk(prev_slns[0], "rho_end3d.vtk", "density", true);


		//Cleanup
			delete mass_matrix;
			delete matrix;			
delete dS_matrix;
delete dg_matrix;
delete K_matrix;





			delete num_flux;






  // Wait for the view to be closed.
  View::wait();
  return 0;
}

