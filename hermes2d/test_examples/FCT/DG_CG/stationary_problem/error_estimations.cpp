/*class Error_WF : public WeakForm<double>
{
public:



  Error_WF(Solution<double>* sln_1,Solution<double>* sln_2 , bool DG)
  : WeakForm<double>(1)
{
 this->set_ext(Hermes::vector<MeshFunction<double>*>(sln_1, sln_2));
 if(DG)
 	add_vector_form_DG(new CustomVectorFormDG);
 else
 	add_vector_form_surf(new Error_WF::CustomVectorFormSurface);

}
  WeakForm<double>* clone() const {return new Error_WF(*this);}


  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface() : VectorFormSurf<double>(0) {};

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
	{
	  double result = 0;
	  for (int i = 0; i < n; i++) 
	  {
		 double v_x = (0.5- e->y[i]); double v_y =(e->x[i]-0.5);
		 double b_dot_n = static_cast<Error_WF*>(wf)->calculate_b_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
		 result += wt[i] *(ext[0]->val[i]-ext[1]->val[i])*(ext[0]->val[i]-ext[1]->val[i])*b_dot_n;
	  }
	  return result;

	}

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
    {
  			return Ord(10);
  	}

    VectorFormSurf<double>* clone() const {return new Error_WF::CustomVectorFormSurface;}


  };
  
    class CustomVectorFormDG: public VectorFormDG<double>
  {
  public:
    CustomVectorFormDG() : VectorFormDG<double>(0) {};

   double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
    {
      double result = 0;
      for (int i = 0;i < n;i++)
      {
       double v_x = (0.5- e->y[i]); double v_y =(e->x[i]-0.5);
		 double b_dot_n = static_cast<Error_WF*>(wf)->calculate_b_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
        //result += wt[i]*((ext[0]->get_val_central(i)-ext[1]->get_val_central(i)) - (ext[0]->get_val_neighbor(i)-ext[1]->get_val_neighbor(i))) * ((ext[0]->get_val_central(i)-ext[1]->get_val_central(i)) - (ext[0]->get_val_neighbor(i)-ext[1]->get_val_neighbor(i)))*b_dot_n;
result += wt[i]*(ext[0]->get_val_central(i) - ext[0]->get_val_neighbor(i)) * (ext[0]->get_val_central(i) - ext[0]->get_val_neighbor(i))*b_dot_n;

//Jump von exakter Loesung = 0, deswegen nur ext[0]!
      }
      return result;

    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
        {
      return Ord(10);
    }

    VectorFormDG<double>* clone() const { return new Error_WF::CustomVectorFormDG; }

  };
  
  
  double calculate_b_dot_v(double x, double y, double vx, double vy) const
  {
		return std::abs(x*vx + y*vy);

	}

  Ord calculate_b_dot_v(Ord x, Ord y, Ord vx, Ord vy) const
  {  
  return Ord(10);
  }

};
*/


  class InterfaceErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
  {
  public:
    InterfaceErrorForm ()  : KellyTypeAdapt<double>::ErrorEstimatorForm(0) 
    {
    this->setAsInterface();
    };

    double value(int n, double *wt, 
               Func<double> *u_ext[], Func<double> *u, 
               Geom<double> *e, Func<double>* *ext) const
	{
	  double result = 0;
	  for (int i = 0; i < n; i++) 
	  {
		 double v_x = (0.5- e->y[i]); double v_y =(e->x[i]-0.5);
		 double b_dot_n = calculate_b_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
result += wt[i]*(u->get_val_central(i) - u->get_val_neighbor(i)) * (u->get_val_central(i) - u->get_val_neighbor(i))*b_dot_n;
	  }
	  return result;

	}

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
    {
  			return Ord(10);
  	}

    
      double calculate_b_dot_v(double x, double y, double vx, double vy) const
  {
		return std::abs(x*vx + y*vy);

	}


  };



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

double calc_error_bdry(Solution<double>* u_1, Solution<double>* u_2,Space<double>* space)
{
 /* Hermes::Mixins::Loggable::Static::info(" error_bdry"); 
	//Error_WF  err_bdry(u_1, u_2, false);
	Error_WF  err_dg(u_1, u_2, true);
	int ref_ndof = space->get_num_dofs();	
	//UMFPackVector<double> rhs_1(ref_ndof); 
	UMFPackVector<double> rhs_2(ref_ndof); 
	//DiscreteProblem<double> dp_1(&err_bdry, space);
	DiscreteProblem<double> dp_2(&err_dg, space);
	//dp_1.assemble(&rhs_1);
	dp_2.assemble(&rhs_2);
double total =0.;
	for(int i = 0; i< (ref_ndof); i++) 
	{
		total += rhs_2.get(i);	
	}
	
		Hermes::Mixins::Loggable::Static::info(" total_error = %f", total);
		return total;*/
		
		
	KellyTypeAdapt<double>* error_estimation = new KellyTypeAdapt<double>(space);
	 error_estimation->add_error_estimator_surf(new InterfaceErrorForm);	
	
double err_est = error_estimation->calc_err_est(u_1);

		Hermes::Mixins::Loggable::Static::info(" total_error = %f", err_est);

		delete error_estimation;
		
		return err_est;
}



double calc_error_streamline_sqr(Solution<double>* u_1, Solution<double>* u_2,Space<double>* space)
{

CustomV_x beta_x(space->get_mesh());
CustomV_y beta_y(space->get_mesh());

      // order of integral
      const int order = 10;

      // initialize total_error
    double total_error = 0.0;
    double beta_K, error_elem;

      // element for the element loop
      Element *e =NULL;

      // refmap for computing Jacobian
			RefMap* rm = new RefMap;
      	rm->set_quad_2d(&g_quad_2d_std);

      // loop over elements
      for_all_active_elements( e, space->get_mesh() ) 
      {
				beta_K =0.;
				error_elem = 0.;
		     // set up the solution quadrature
		     u_1->set_quad_2d(&g_quad_2d_std);
		     u_1->set_active_element(e);
		     u_1->set_quad_order(order);		     
		     u_2->set_quad_2d(&g_quad_2d_std);
		     u_2->set_active_element(e);
		     u_2->set_quad_order(order);	
		     beta_x.set_quad_2d(&g_quad_2d_std);
		     beta_x.set_active_element(e);
		     beta_x.set_quad_order(order);
		     beta_y.set_quad_2d(&g_quad_2d_std);
		     beta_y.set_active_element(e);
		     beta_y.set_quad_order(order);			

        // get the constant Jacobian
        rm->set_active_element(e);
        double jac = rm->get_const_jacobian();

        // get the quadrature points
        int np = u_1->get_quad_2d()->get_num_points(order,HERMES_MODE_QUAD);
        double3 *pt = g_quad_2d_std.get_points(order,HERMES_MODE_QUAD);

		     // get the function derivative values
		     Func<double>* u = init_fn( u_1, order );
			  Func<double>* v = init_fn( u_2, order );
			  Func<double>* beta_1 = init_fn( &beta_x, order );
			  Func<double>* beta_2 = init_fn( &beta_y, order );

        // loop over points and integrate the energy
        for( int j = 0; j < np; ++j ) {
					error_elem +=pt[j][2]*jac* ((u->dx[j]-v->dx[j])*beta_1->val[j] + (u->dy[j]-v->dy[j])*beta_2->val[j] );   
					beta_K += pt[j][2]*jac*(beta_1->val[j]*beta_1->val[j] + beta_2->val[j]*beta_2->val[j]);
        }
        		
				total_error += error_elem* e->get_diameter()/std::sqrt(beta_K);
				v->free_fn();
				u->free_fn();
				beta_1->free_fn();
				beta_2->free_fn();
				delete beta_1;
				delete beta_2;
				delete v;
				delete u;
      }
		delete rm;

}
