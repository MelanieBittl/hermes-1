
void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
int ndof = space->get_num_dofs();

  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  errorCalculator_l2.add_error_form(new CustomNormFormVol(0,0));
  errorCalculator_l2.calculate_errors(u_new, u_prev_time,false);

double err_l2_2 = errorCalculator_l2.get_total_error_squared();


Hermes::Mixins::Loggable::Static::info("l2=%.3e,  ndof = %d",Hermes::sqrt(err_l2_2),  ndof);

FILE * pFile;
pFile = fopen ("l2.dat","a");
    fprintf (pFile, "%d	%.4e \n", ndof, Hermes::sqrt(err_l2_2));
fclose (pFile);


// Output solution in VTK format.
/*	Linearizer lin;
MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  */




}