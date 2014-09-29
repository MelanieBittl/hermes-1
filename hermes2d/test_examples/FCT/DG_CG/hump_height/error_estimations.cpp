
void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space, bool end = true)
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


// Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
if(end)
{
lin.save_solution_vtk(u_new, "sln.vtk", "solution", mode_3D);
//lin.save_solution_vtk(u_prev_time, "exact.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  
}else{

lin.save_solution_vtk(u_new, "sln_mid.vtk", "solution", mode_3D);
lin.save_solution_vtk(u_prev_time, "exact_mid.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error_mid.vtk" , "error", false);  
}


}
