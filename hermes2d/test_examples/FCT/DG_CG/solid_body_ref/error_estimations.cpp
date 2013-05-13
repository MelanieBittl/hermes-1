

double calc_error_l2(Solution<double>* u_1, Solution<double>* u_2,Space<double>* space)
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
		     // set up the solution quadrature
		     u_1->set_quad_2d(&g_quad_2d_std);
		     u_1->set_active_element(e);
		     u_1->set_quad_order(order);		     
		     u_2->set_quad_2d(&g_quad_2d_std);
		     u_2->set_active_element(e);
		     u_2->set_quad_order(order);		    

		     // get the quadrature points
		     int np = u_1->get_quad_2d()->get_num_points(order,HERMES_MODE_QUAD);
		     double3 *pt = g_quad_2d_std.get_points(order,HERMES_MODE_QUAD);
		             // get the constant Jacobian
        rm->set_active_element(e);
        double jac = rm->get_const_jacobian();

		     // get the function derivative values
		     Func<double>* u = init_fn( u_1, order );
			  Func<double>* v = init_fn( u_2, order );

		     // loop over points and determine err_max
		     for( int j = 0; j < np; ++j ) 
					err_total +=pt[j][2]*jac*Hermes::sqr(u->val[j] - v->val[j]) ;          
		     

				v->free_fn();
				u->free_fn();
				delete v;
				delete u;
		  }
		  delete rm;
	return std::sqrt(err_total);

}



double calc_error_max(Solution<double>* u_1, Solution<double>* u_2,Space<double>* space)
{

		// order of integral
		const int order = 10;
		double err_max =0.;

		// element for the element loop
		Element *e =NULL;
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
		     int np = u_1->get_quad_2d()->get_num_points(order,HERMES_MODE_QUAD);
		     double3 *pt = g_quad_2d_std.get_points(order,HERMES_MODE_QUAD);

		     // get the function derivative values
		     Func<double>* u = init_fn( u_1, order );
			  Func<double>* v = init_fn( u_2, order );

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

