
void calc_error_total_new(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
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
pFile = fopen ("error.txt","a");
    fprintf (pFile, "l2=%.4e, h1=%.4e, h1_semi= %.4e,ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_h1_2),Hermes::sqrt(err_h1_semi_2), ndof);
fclose (pFile); 




// Output solution in VTK format.
	Linearizer lin;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  
}
