#define HERMES_REPORT_ALL
#include "definitions.h"
#include "hermes2d.h"
using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT =2;       						// Initial polynomial degree.



MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";



double calc_error_l2(MeshFunctionSharedPtr<double> u_1, MeshFunctionSharedPtr<double> u_2,SpaceSharedPtr<double> space)
{
Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_TRIANGLE;
//Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;
	// order of integral
	const int order = 10;
	double err_total =0.;
	// element for the element loop
	Element *e =NULL;
	// refmap for computing Jacobian
	RefMap* rm = new RefMap;
	rm->set_quad_2d(&g_quad_2d_std);

	for_all_active_elements(e, space->get_mesh())
	{
			// set up the solution quadrature
			u_1->set_quad_2d(&g_quad_2d_std);
			u_1->set_active_element(e);
			u_1->set_quad_order(order);	
			u_2->set_quad_2d(&g_quad_2d_std);
			u_2->set_active_element(e);
			u_2->set_quad_order(order);	

			// get the quadrature points
			int np = u_1->get_quad_2d()->get_num_points(order,mode);
			double3 *pt = g_quad_2d_std.get_points(order, mode);
			// get the constant Jacobian
		  rm->set_active_element(e);
		  double jac = rm->get_const_jacobian();
			MeshFunction<double>* sln = u_1->clone();
			sln->set_active_element(e);
			MeshFunction<double>* ref_sln = u_2->clone();
			ref_sln->set_active_element(e);

			// get the function derivative values
			Func<double>* u = init_fn( sln, order );
			Func<double>* v = init_fn( ref_sln, order );
		double diam = e->get_diameter();
		double area =Hermes::sqrt(e->get_area());
		double v_x = 0.5; 
		double v_y = 1.;
		double abs_v = Hermes::sqrt(v_x*v_x+v_y*v_y);

		
			for( int j = 0; j < np; ++j )
				//err_total +=pt[j][2]*jac*Hermes::sqr(v_x*(u->dx[j])+v_y*(u->dy[j]))*diam/(abs_v*area);
				err_total += pt[j][2]*jac*Hermes::sqr(v_x*(u->dx[j]-v->dx[j])+v_y*(u->dy[j]-v->dy[j]))*diam/(abs_v*area);

			v->free_fn();
			u->free_fn();
			delete v;
			delete u;
			delete sln;
			delete ref_sln;
	}
	delete rm;
	return Hermes::sqrt(err_total);
}



int main(int argc, char* argv[])
{ 
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("unit.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements(); 

//Hermes2DApi.set_integral_param_value(numThreads, 1); 

  // Create an L2 	space with default shapeset.
//CustomDirichletCondition bc_essential(Hermes::vector<std::string>("inlet1","inlet2"));
 // EssentialBCs<double>  bcs(&bc_essential);
  SpaceSharedPtr<double> space(new L2_SEMI_CG_Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new L2Space<double>(mesh,P_INIT));	
  //SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));	

/*BaseView<double> bview("Baseview", new WinGeom(450, 0, 440, 350));
bview.show(space);
View::wait(HERMES_WAIT_KEYPRESS);*/

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>  u_new(new Solution<double>);
  MeshFunctionSharedPtr<double> u_prev_time(new CustomInitialCondition(mesh));
 
  // Initialize views.
	ScalarView sview("solution_1", new WinGeom(500, 500, 500, 400));
	ScalarView s2view("solution_2", new WinGeom(0, 500, 500, 400));

	ScalarView fview("filter", new WinGeom(500, 500, 500, 400));
	ScalarView lview("initial condition", new WinGeom(500, 0, 500, 400));
	lview.show(u_prev_time);
fview.set_min_max_range(0,0.001);
//View::wait(HERMES_WAIT_KEYPRESS);
	
	int ndof = space->get_num_dofs();


	double *coeff_vec = new double[ndof];		
	double *coeff_vec_2 = new double[ndof];
	double *coeff_vec_3 = new double[ndof];
	double *coeff_vec_4 = new double[ndof];
	double* vec_new;
  memset(coeff_vec_2, 0, ndof*sizeof(double));
 memset(coeff_vec_3, 0, ndof*sizeof(double));
 memset(coeff_vec_4, 0, ndof*sizeof(double));


	///////////////////////////------------false, false (only CG), false, true (DG) --------------------------------------------------------------------------------
	CustomWeakForm wf_surf(u_prev_time,true, true);
	UMFPackMatrix<double>* dg_surface_matrix = new UMFPackMatrix<double> ; 
	UMFPackVector<double> * surf_rhs = new UMFPackVector<double> (ndof); 
	DiscreteProblem<double> * dp_surf = new DiscreteProblem<double> (&wf_surf,space);	
	dp_surf->set_linear(true,false);
	dp_surf->assemble(dg_surface_matrix,surf_rhs);
 UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>( dg_surface_matrix, surf_rhs);    
     if(solver->solve())
     {
     	vec_new = solver->get_sln_vector();
      Solution<double>::vector_to_solution(vec_new, space, u_new);
	for(int i=0; i<ndof; i++) coeff_vec_2[i] = vec_new[i];
      }
    else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		sview.show(u_new);


/*
CustomWeakForm wf(u_prev_time,true, false);
 // Initialize linear solver.
 Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();
	for(int i=0; i<ndof; i++) coeff_vec_2[i] = sln_vector[i];

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, u_new);
		s2view.show(u_new);
  }catch(std::exception& e)
  {
    std::cout << e.what();
  }
*/

 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(mesh));
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG(AbsoluteError);
  ErrorCalculator<double> errorCalculator_sd(AbsoluteError);
  errorCalculator_l2.add_error_form(new CustomNormFormVol(0,0));
	errorCalculator_sd.add_error_form(new StreamlineDiffusionNorm(0,0));
  errorCalculator_surf.add_error_form(new CustomNormFormSurf(0,0));
  errorCalculator_DG.add_error_form(new CustomNormFormDG(0,0));

  errorCalculator_l2.calculate_errors(u_new, u_prev_time);
	errorCalculator_surf.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG.calculate_errors(u_new,u_prev_time);
	errorCalculator_sd.calculate_errors(u_new, u_prev_time);


double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_surf_2 = errorCalculator_surf.get_total_error_squared();
double err_DG_2 = errorCalculator_DG.get_total_error_squared();
double err_sd_2 = errorCalculator_sd.get_total_error_squared();




double test = calc_error_l2(u_new, u_prev_time,space);
Hermes::Mixins::Loggable::Static::info("test=%.3e", test);

double total = Hermes::sqrt(err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2);

Hermes::Mixins::Loggable::Static::info("l2=%.3e, surf = %.3e, dg = %.3e, sd = %.3e, total= %.3e, ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof);

//MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));

//fview.show(filter);
/*
  OGProjection<double> ogProjection;	
        ogProjection.project_global(space, u_prev_time, coeff_vec_3, HERMES_L2_NORM);	

 Wf_residual wf_test(u_new, u_prev_time);
	UMFPackVector<double> * rhs = new UMFPackVector<double> (ndof); 
UMFPackMatrix<double>* matrix = new UMFPackMatrix<double> ; 
	DiscreteProblem<double> * dp_test = new DiscreteProblem<double> (&wf_test,space);	
	dp_test->set_linear(true,false);
	dp_test->assemble(matrix,rhs);

dg_surface_matrix->multiply_with_vector(coeff_vec_3, coeff_vec_4);
double test =0; double test_2 =0;

for(int i =0; i<ndof;i++) 
{
	test+=rhs->get(i);
test_2 += (coeff_vec_4[i]-surf_rhs->get(i));
}

Hermes::Mixins::Loggable::Static::info("l(u)=%.5e, A_u_init-f=%.5e",test, test_2);
	*/




  // Wait for the view to be closed.
  View::wait();
  return 0;
}

