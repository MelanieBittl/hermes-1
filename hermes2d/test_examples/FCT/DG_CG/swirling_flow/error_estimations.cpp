void calc_error_total(MeshFunctionSharedPtr<double> u_new,MeshFunctionSharedPtr<double> u_proj, MeshFunctionSharedPtr<double> u_prev_time, SpaceSharedPtr<double> space, double final_time, double time)
{
int ndof = space->get_num_dofs();
 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(space->get_mesh()));
 
  ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG(AbsoluteError);
  ErrorCalculator<double> errorCalculator_sd(AbsoluteError);
  errorCalculator_l2.add_error_form(new CustomNormFormVol(0,0));
	errorCalculator_sd.add_error_form(new StreamlineDiffusionNorm(0,0,space->get_mesh(),final_time,time));
  errorCalculator_surf.add_error_form(new CustomNormFormSurf(0,0,final_time,time));
  errorCalculator_DG.add_error_form(new CustomNormFormDG(0,0,final_time,time));

  errorCalculator_l2.calculate_errors(u_new, u_prev_time);
	errorCalculator_surf.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG.calculate_errors(u_new,u_prev_time);
	errorCalculator_sd.calculate_errors(u_new, u_prev_time);


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



double total = Hermes::sqrt(err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2);

Hermes::Mixins::Loggable::Static::info("init: l2=%.3e, surf = %.3e, dg = %.3e, sd = %.3e, total= %.3e, ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof);
/*
FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%.4e, surf = %.4e, dg = %.4e, sd = %.4e, total= %.4e, ndof = %d, diam_min =%.4e, diam_max =%.4e",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total, ndof,diam_min, diam_max);
fclose (pFile); */

FILE * pFile_2;
pFile_2 = fopen ("error_s_init.txt","a");
    fprintf (pFile_2, "l2=%.4e, total= %.4e, ndof = %d, diam_ =%.4e \n",
Hermes::sqrt(err_l2_2), total , ndof,diam_max);
fclose (pFile_2); 




//for projected_method
  errorCalculator_l2.calculate_errors(u_new, u_proj);
	errorCalculator_surf.calculate_errors(u_new, u_proj);	
		errorCalculator_DG.calculate_errors(u_new,u_proj);
	errorCalculator_sd.calculate_errors(u_new, u_proj);


 err_l2_2 = errorCalculator_l2.get_total_error_squared();
 err_surf_2 = errorCalculator_surf.get_total_error_squared();
 err_sd_2 = errorCalculator_sd.get_total_error_squared();
 err_DG_2 = errorCalculator_DG.get_total_error_squared();

total = Hermes::sqrt(err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2);

Hermes::Mixins::Loggable::Static::info("proj: l2=%.3e, surf = %.3e, dg = %.3e, sd = %.3e, total= %.3e, ndof = %d",
Hermes::sqrt(err_l2_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof);

pFile_2 = fopen ("error_s_proj.txt","a");
    fprintf (pFile_2, "l2=%.4e, total= %.4e, ndof = %d, diam_ =%.4e \n",
Hermes::sqrt(err_l2_2), total , ndof,diam_max);
fclose (pFile_2); 

// Output solution in VTK format.
Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
lin.save_solution_vtk(u_new, "sln.vtk", "solution", mode_3D);
//lin.save_solution_vtk(u_prev_time, "init.vtk", "solution", mode_3D);
lin.save_solution_vtk(filter, "error.vtk" , "error", false);  




}
