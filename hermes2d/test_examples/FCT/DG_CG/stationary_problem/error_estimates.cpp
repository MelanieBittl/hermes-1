
double calc_error_l2(MeshFunctionSharedPtr<double> u_1, MeshFunctionSharedPtr<double> u_2,SpaceSharedPtr<double> space)
{
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
Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;
			if(e->is_triangle()) mode = HERMES_MODE_TRIANGLE;
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

		double err_elem = 0.;	
		double abs_v =0;	

		double* x_coord = rm->get_phys_x(order);
		double* y_coord = rm->get_phys_y(order);
			for( int j = 0; j < np; ++j )
			{
				double v_x = 0.5;
				double v_y = 1.;
 				//double v_x = y_coord[j]; 
				//double v_y = 1.-x_coord[j];
				err_elem += pt[j][2]*jac*Hermes::sqr(v_x*(u->dx[j]-v->dx[j])+v_y*(u->dy[j]-v->dy[j]));
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
//err_elem += pt[j][2]*jac*Hermes::sqr(u->val[j]-v->val[j]);
			}
			err_total += (err_elem*diam/Hermes::sqrt(abs_v));
			//err_total+= err_elem;
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

void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
int ndof = space->get_num_dofs();
 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(space->get_mesh()));
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG(AbsoluteError);
  ErrorCalculator<double> errorCalculator_sd(AbsoluteError);
  errorCalculator_l2.add_error_form(new CustomNormFormVol(0,0));
	errorCalculator_sd.add_error_form(new StreamlineDiffusionNorm(0,0,space->get_mesh()));
  errorCalculator_surf.add_error_form(new CustomNormFormSurf(0,0));
  errorCalculator_DG.add_error_form(new CustomNormFormDG(0,0));

  errorCalculator_l2.calculate_errors(u_new, u_prev_time,false);
	errorCalculator_surf.calculate_errors(u_new, u_prev_time,false);
	errorCalculator_DG.calculate_errors(u_new,sln_zero, false);
	errorCalculator_sd.calculate_errors(u_new, u_prev_time, false);


double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_surf_2 = errorCalculator_surf.get_total_error_squared();
double err_DG_2 = errorCalculator_DG.get_total_error_squared();
double err_sd_2 = errorCalculator_sd.get_total_error_squared();

double diam;Element* e;
double diam_max =0.; double diam_min = 100.;
for_all_active_elements(e, space->get_mesh()) 
{	
	diam = e->get_diameter(); // Laenge der kleinste Kante
	if(diam<diam_min)	diam_min = diam;
	if(diam>diam_max) diam_max = diam;
	break;
}

//double test = calc_error_l2(u_new, u_prev_time,space);
//Hermes::Mixins::Loggable::Static::info("test=%.3e", test);

double total = Hermes::sqrt(err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2);

Hermes::Mixins::Loggable::Static::info("l2=%.3e, surf = %.3e, dg = %.3e, sd = %.3e, total= %.3e, ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof);

FILE * pFile;
pFile = fopen ("error.txt","a");
    fprintf (pFile, "l2=%.4e, surf = %.4e, dg = %.4e, sd = %.4e, total= %.4e, ndof = %d, diam_min =%.4e, diam_max =%.4e \n",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof,diam_min, diam_max);
fclose (pFile);  

FILE * pFile_2;
pFile_2 = fopen ("error_s.txt","a");
    fprintf (pFile_2, "l2=%.4e, total= %.4e, ndof = %d, diam_ =%.4e \n",
Hermes::sqrt(err_l2_2), total , ndof,diam_max);
fclose (pFile_2); 


pFile_2 = fopen ("error_l2.txt","a");
    fprintf (pFile_2, "	%.4e \n", Hermes::sqrt(err_l2_2));
fclose (pFile_2); 

pFile_2 = fopen ("error_total.txt","a");
    fprintf (pFile_2, "	%.4e \n", total);
fclose (pFile_2); 


// Output solution in VTK format.

Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//ScalarView fview("filter", new WinGeom(500, 500, 500, 400));
//fview.show(filter);
lin.save_solution_vtk(u_new, "sln.vtk", "solution", mode_3D);
//lin.save_solution_vtk(u_prev_time, "init.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  




}
