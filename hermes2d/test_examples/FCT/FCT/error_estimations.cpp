
 
double* calc_error( SpaceSharedPtr<double> ref_space, MeshFunctionSharedPtr<double> u_fct,MeshFunctionSharedPtr<double> u_init)
{
      // order of integral
      const int order = 10;

      // initialize total_error
      double* total_error = new double[2];
			total_error[0] = 0.;
			total_error[1] = 0.;

      // element for the element loop
      Element *e =NULL;

      // refmap for computing Jacobian
			RefMap* rm = new RefMap;
      rm->set_quad_2d(&g_quad_2d_std);

Solution<double>* u = new Solution<double>; u = static_cast<Solution<double>* > (u_fct->clone());
Solution<double>* u_exact = new Solution<double>; u_exact = static_cast<Solution<double>* > (u_init->clone());


      // loop over elements
      for_all_active_elements( e, ref_space->get_mesh() ) {


        // set up the solution quadrature
        u->set_quad_2d(&g_quad_2d_std);
        u->set_active_element(e);
        u->set_quad_order(order);
				
				u_exact->set_quad_2d(&g_quad_2d_std);
        u_exact->set_active_element(e);
       u_exact->set_quad_order(order);


        // get the constant Jacobian
        rm->set_active_element(e);
        double jac = rm->get_const_jacobian();

        // get the quadrature points
        int np = u->get_quad_2d()->get_num_points(order,HERMES_MODE_QUAD);
        double3 *pt = g_quad_2d_std.get_points(order,HERMES_MODE_QUAD);

        // get the function derivative values
        Func<double>* v = init_fn( u, order );
				Func<double>* v_ex = init_fn( u_exact, order );

		double* x_coord = rm->get_phys_x(order);
		double* y_coord = rm->get_phys_y(order);

        // loop over points and integrate the energy
        for( int j = 0; j < np; ++j ) {
					if((x_coord[j]>=0.375)&&(x_coord[j]<=0.625))
					{
							if((y_coord[j]>=0.375)&&(y_coord[j]<=0.625))
									total_error[0] +=pt[j][2]*jac*( (v->val[j] - v_ex->val[j])*(v->val[j] - v_ex->val[j]) );   
					} 
					if((x_coord[j]>=0.4375)&&(x_coord[j]<=0.5625))
					{
							if((y_coord[j]>=0.4375)&&(y_coord[j]<=0.5625))
									total_error[1] +=pt[j][2]*jac*( (v->val[j] - v_ex->val[j])*(v->val[j] - v_ex->val[j]) );   
					}      
        }


			v->free_fn();
			v_ex->free_fn();

				delete v;
			delete v_ex;
      }
		delete rm;
	delete u; 
	delete u_exact;

return total_error;
 }


void calc_error_total_time(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space, double time)
{
int ndof = space->get_num_dofs();
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
  errorCalculator_l2.calculate_errors(u_new, u_prev_time,false);
double err_l2_2 = errorCalculator_l2.get_total_error_squared();
//double* err_limit = calc_error(space, u_new, u_prev_time);


FILE * pFile;
pFile = fopen ("error.txt","a");
fprintf (pFile, "time =%g:  l2=%.4e, ndof = %d \n",time, Hermes::sqrt(err_l2_2),ndof);
    //fprintf (pFile, "time =%g:  l2=%.4e, limit_1=%.4e, limit_2=%.4e, ndof = %d",time, Hermes::sqrt(err_l2_2),Hermes::sqrt(err_limit[0]),Hermes::sqrt(err_limit[1]),ndof);
fclose (pFile); 

//delete [] err_limit;

}

void calc_error_proj(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
	int ndof = space->get_num_dofs();
	 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_h1(AbsoluteError);
    ErrorCalculator<double> errorCalculator_h1_semi(AbsoluteError);

  errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
    errorCalculator_h1.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_NORM));
  errorCalculator_h1_semi.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_SEMINORM));

  errorCalculator_l2.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1_semi.calculate_errors(u_new, u_prev_time);


double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_h1_2 = errorCalculator_h1.get_total_error_squared();
double err_h1_semi_2 = errorCalculator_h1.get_total_error_squared();


Hermes::Mixins::Loggable::Static::info("l2=%.4e, ndof = %d", Hermes::sqrt(err_l2_2),ndof);

FILE * pFile;
pFile = fopen ("error_proj.txt","w");
    fprintf (pFile, "l2=%.4e, h1=%.4e, h1_semi= %.4e,ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_h1_2),Hermes::sqrt(err_h1_semi_2), ndof);
fclose (pFile); 

}


void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
int ndof = space->get_num_dofs();
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_h1(AbsoluteError);
    ErrorCalculator<double> errorCalculator_h1_semi(AbsoluteError);

  errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
    errorCalculator_h1.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_NORM));
  errorCalculator_h1_semi.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_SEMINORM));

  errorCalculator_l2.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1.calculate_errors(u_new, u_prev_time);
  errorCalculator_h1_semi.calculate_errors(u_new, u_prev_time);


double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_h1_2 = errorCalculator_h1.get_total_error_squared();
double err_h1_semi_2 = errorCalculator_h1.get_total_error_squared();


Hermes::Mixins::Loggable::Static::info("l2=%.4e, ndof = %d", Hermes::sqrt(err_l2_2),ndof);

FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%.4e, h1=%.4e, h1_semi= %.4e,ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_h1_2),Hermes::sqrt(err_h1_semi_2), ndof);
fclose (pFile); 




// Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
lin.save_solution_vtk(u_new, "sln3d.vtk", "solution", mode_3D);
lin.save_solution_vtk(u_new, "sln2d.vtk", "solution", false);
//lin.save_solution_vtk(u_prev_time, "init.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  

}



