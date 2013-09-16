double calc_error_max(MeshFunctionSharedPtr<double> u_1, MeshFunctionSharedPtr<double> u_2,SpaceSharedPtr<double> space)
{

Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;

		// order of integral
		const int order = 10;
		double err_max =0.;

		// element for the element loop
		Element *e =NULL;
		for_all_active_elements(e, space->get_mesh())
		{
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
		     double3 *pt = g_quad_2d_std.get_points(order,mode);

			MeshFunction<double>* sln = u_1->clone();
			sln->set_active_element(e);
			MeshFunction<double>* ref_sln = u_2->clone();
			ref_sln->set_active_element(e);

			// get the function derivative values
			Func<double>* u = init_fn( sln, order );
			Func<double>* v = init_fn( ref_sln, order );

		     // loop over points and determine err_max
		     for( int j = 0; j < np; ++j ) 
					if(err_max<	std::abs(u->val[j] - v->val[j]))
						err_max = std::abs(u->val[j] - v->val[j]);         
		     

				v->free_fn();
				u->free_fn();
				delete v;
				delete u;
		  }
	return err_max;
}

double calc_error_max(double* u, double* v, int ndof)
{
	double max =0;	
	for(int i =0; i<ndof; i++)
	{
		if(max< std::abs(u[i]-v[i]))
			max = std::abs(u[i]-v[i]);	
	}

	return max;
}


double calc_coord( SpaceSharedPtr<double> ref_space, MeshFunctionSharedPtr<double> u_fct, int comp)
{
      // order of integral
      const int order = 10;

			double result = 0.;
      // element for the element loop
      Element *e =NULL;

      // refmap for computing Jacobian
			RefMap* rm = new RefMap;
      rm->set_quad_2d(&g_quad_2d_std);

			Solution<double>* u = new Solution<double>;
				u =  static_cast<Solution<double>* > (u_fct->clone());

			Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;
      // loop over elements
      for_all_active_elements( e, ref_space->get_mesh() ) 
			{				
				if(e->is_triangle()) mode = HERMES_MODE_TRIANGLE;
        // set up the solution quadrature
        u->set_quad_2d(&g_quad_2d_std);
        u->set_active_element(e);
        u->set_quad_order(order);

        // get the constant Jacobian
        rm->set_active_element(e);
        double jac = rm->get_const_jacobian();

        // get the quadrature points
        int np = u->get_quad_2d()->get_num_points(order,mode);
        double3 *pt = g_quad_2d_std.get_points(order,mode);

				double* x =  rm->get_phys_x(order);
				double* y =  rm->get_phys_y(order);

        // get the function derivative values
        Func<double>* v = init_fn( u, order );

        // loop over points and integrate the energy
        for( int j = 0; j < np; ++j ) 
				{
					if(comp==0)
						result +=pt[j][2]*jac*(x[j]*v->val[j]);  
					if(comp==1)      
						result +=pt[j][2]*jac*(y[j]*v->val[j]);  
        }

				v->free_fn();
				delete v;

      }
		delete rm;
		delete u; 

	return result;
 }





void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
int ndof = space->get_num_dofs();
 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(space->get_mesh()));
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_h1(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG(AbsoluteError);
  ErrorCalculator<double> errorCalculator_sd(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf_1(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG_1(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf_2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG_2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_h1_semi(AbsoluteError);

  errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
  errorCalculator_h1.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_NORM));
  errorCalculator_h1_semi.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_SEMINORM));

	errorCalculator_sd.add_error_form(new StreamlineDiffusionNorm(0,0,space->get_mesh()));
  errorCalculator_surf.add_error_form(new CustomNormFormSurf(0,0));
  errorCalculator_DG.add_error_form(new CustomNormFormDG(0,0));
  errorCalculator_surf_1.add_error_form(new CustomNormFormSurf_1(0,0));
  errorCalculator_DG_1.add_error_form(new CustomNormFormDG_1(0,0));
  errorCalculator_surf_2.add_error_form(new CustomNormFormSurf_2(0,0));
  errorCalculator_DG_2.add_error_form(new CustomNormFormDG_2(0,0));


  errorCalculator_l2.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1_semi.calculate_errors(u_new, u_prev_time);
	errorCalculator_sd.calculate_errors(u_new, u_prev_time);
	errorCalculator_surf.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG.calculate_errors(u_new,sln_zero);
	errorCalculator_surf_1.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG_1.calculate_errors(u_new,sln_zero);
	errorCalculator_surf_2.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG_2.calculate_errors(u_new,u_prev_time);



double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_h1_2 = errorCalculator_h1.get_total_error_squared();
double err_h1_semi_2 = errorCalculator_h1.get_total_error_squared();
double err_surf_2 = errorCalculator_surf.get_total_error_squared();
double err_DG_2 = errorCalculator_DG.get_total_error_squared();
double err_sd_2 = errorCalculator_sd.get_total_error_squared();
double err_surf_1_2 = errorCalculator_surf_1.get_total_error_squared();
double err_DG_1_2 = errorCalculator_DG_1.get_total_error_squared();
double err_surf_2_2 = errorCalculator_surf_2.get_total_error_squared();
double err_DG_2_2 = errorCalculator_DG_2.get_total_error_squared();

double diam;Element* e;
double diam_max =0.; double diam_min = 100.;
for_all_active_elements(e, space->get_mesh()) 
{	
	diam = e->get_diameter(); // Laenge der kleinste Kante
	if(diam<diam_min)	diam_min = diam;
	if(diam>diam_max) diam_max = diam;
	break;
}

/*
double x_h = calc_coord(space, u_new, 0);
double y_h = calc_coord(space, u_new, 1);

Hermes::Mixins::Loggable::Static::info("x_h=%f, y_h = %f,", x_h, y_h);*/


//double err_max = calc_error_max(u_new, u_prev_time,space);

double total_conv = err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2;
double total_diff = err_h1_semi_2+err_surf_1_2+err_surf_2_2+err_DG_1_2+err_DG_2_2;

double total = Hermes::sqrt(total_conv+total_diff);

Hermes::Mixins::Loggable::Static::info("l2=%.3e, h1=%.3e, total= %.3e, ndof = %d", 
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_h1_2), total , ndof);

FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%.4e, h1=%.4e, diff= %.4e, conv = %.4e, total= %.4e, ndof = %d, diam =%.4e",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_h1_2),Hermes::sqrt(total_diff),Hermes::sqrt(total_conv), total, ndof,diam_max);
fclose (pFile);  

/*
// Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
lin.save_solution_vtk(u_new, "sln.vtk", "solution", mode_3D);
lin.save_solution_vtk(u_prev_time, "exact.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  
*/


}
