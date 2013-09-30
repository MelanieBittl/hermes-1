void calculate_osc(MeshFunctionSharedPtr<double> u_new,SpaceSharedPtr<double> space, double* result)
{
	Element* e;
double total_int =0.; 
double total_exp =0.;

	for_all_active_elements(e, space->get_mesh()) 
	{
			Node* vn; double xi1 = 0; double xi2 = 0;
			for(int i =0; i<e->get_nvert(); i++)
			{		vn = e->vn[i];
				if((vn->x <=0.5) && (vn->y >= 0.1)) 
					{ 
						u_new->get_refmap()->untransform(e, vn->x, vn->y, xi1, xi2);
						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
						total_int += Hermes::sqr(std::min(func_val,0.)) + Hermes::sqr(std::max (0., func_val-1.));
					}else if(vn->x>=0.7)
					{
						u_new->get_refmap()->untransform(e, vn->x, vn->y, xi1, xi2);
						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
							total_exp += Hermes::sqr(std::max (0., func_val-1.));				
					}
			}
		/*	double x_c, y_c;
			e->get_center(x_c, y_c);
							if((x_c <=0.5) && (y_c >= 0.1)) 
					{ 
						u_new->get_refmap()->untransform(e, x_c, y_c, xi1, xi2);

						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
						total_int += Hermes::sqr(std::min(func_val,0.)) + Hermes::sqr(std::max (0., func_val-1.));
					}else if(x_c >=0.7)
					{
						u_new->get_refmap()->untransform(e, x_c, y_c, xi1, xi2);
						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
							total_exp += Hermes::sqr(std::max (0., func_val-1.));
				
					}*/
	}
	result[0] =Hermes::sqrt(total_int); result[1] =Hermes::sqrt(total_exp);

}

void calculate_smear(MeshFunctionSharedPtr<double> u_new,SpaceSharedPtr<double> space,double* result)
{
	Element* e;
double total_int =0.; 
double total_exp =0.;
double x_1 = 1.; double x_2=1.;


	for_all_active_elements(e, space->get_mesh()) 
	{
			Node* vn; double xi1 = 0; double xi2 = 0;
			for(int i =0; i<e->get_nvert(); i++)
			{		vn = e->vn[i];
				if((vn->y == 0.25)) 
					{ 
						u_new->get_refmap()->untransform(e, vn->x, vn->y, xi1, xi2);
						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
						if((func_val >= 0.1)&&(vn->x < x_1)) x_1 = vn->x;
						if((func_val >= 0.9)&&(vn->x < x_2)) x_2 = vn->x;
					}else if(vn->x>=0.7)
					{
						u_new->get_refmap()->untransform(e, vn->x, vn->y, xi1, xi2);
						double func_val =  (dynamic_cast<Solution<double>*>(u_new.get()))->get_ref_value(e,xi1, xi2);
							total_exp += Hermes::sqr(std::min(func_val-1.,0.));
				
					}
			}



	}


	result[0] = x_2 - x_1; result[1] =Hermes::sqrt(total_exp);


}



void calc_error_total(MeshFunctionSharedPtr<double> u_new, MeshFunctionSharedPtr<double> u_prev_time,SpaceSharedPtr<double> space)
{
int ndof = space->get_num_dofs();
 MeshFunctionSharedPtr<double> sln_zero(new ZeroSolution<double>(space->get_mesh()));


double osc[2]; calculate_osc(u_new,space, osc);
double smear[2]; calculate_smear(u_new,space,smear);


double diam;Element* e;
double diam_max =0.; double diam_min = 100.;
for_all_active_elements(e, space->get_mesh()) 
{	
	diam = e->get_diameter(); // Laenge der kleinste Kante
	if(diam<diam_min)	diam_min = diam;
	if(diam>diam_max) diam_max = diam;
	break;
}

Hermes::Mixins::Loggable::Static::info("osc_int = %.3e, osc_exp = %.3e, smear_int = %.3e, smear_exp = %.3e, diam=%f", osc[0], osc[1], smear[0], smear[1], diam_min); 
FILE * pFile;
pFile = fopen ("osc_smear.txt","w");
    fprintf (pFile, "osc_int = %.3e, osc_exp = %.3e, smear_int = %.3e, smear_exp = %.3e, diam=%f", osc[0], osc[1], smear[0], smear[1], diam_min);
fclose (pFile);  

 /* 
	ErrorCalculator<double> errorCalculator_l2(AbsoluteError);
  ErrorCalculator<double> errorCalculator_h1(AbsoluteError);
  ErrorCalculator<double> errorCalculator_surf(AbsoluteError);
  ErrorCalculator<double> errorCalculator_DG(AbsoluteError);
  ErrorCalculator<double> errorCalculator_sd(AbsoluteError);

  errorCalculator_l2.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_L2_NORM));
  errorCalculator_h1.add_error_form(new DefaultNormFormVol<double>(0,0,HERMES_H1_NORM));
	errorCalculator_sd.add_error_form(new StreamlineDiffusionNorm(0,0,space->get_mesh()));
  errorCalculator_surf.add_error_form(new CustomNormFormSurf(0,0));
  errorCalculator_DG.add_error_form(new CustomNormFormDG(0,0));
  errorCalculator_l2.calculate_errors(u_new, u_prev_time,false);
  errorCalculator_h1.calculate_errors(u_new, u_prev_time,false);
	errorCalculator_surf.calculate_errors(u_new, u_prev_time);
	errorCalculator_DG.calculate_errors(u_new,sln_zero);
	errorCalculator_sd.calculate_errors(u_new, u_prev_time);

double err_l2_2 = errorCalculator_l2.get_total_error_squared();
double err_H1_2 = errorCalculator_h1.get_total_error_squared();
double err_surf_2 = errorCalculator_surf.get_total_error_squared();
double err_DG_2 = errorCalculator_DG.get_total_error_squared();
double err_sd_2 = errorCalculator_sd.get_total_error_squared();

double total = Hermes::sqrt(err_l2_2+0.5*err_surf_2+0.5*err_DG_2+err_sd_2);

Hermes::Mixins::Loggable::Static::info(" l2=%.3e,h1 = %.3e, surf = %.3e, dg = %.3e, sd = %.3e, total= %.3e, ndof = %d",
Hermes::sqrt(err_l2_2),Hermes::sqrt(err_H1_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total , ndof);

//FILE * pFile;
pFile = fopen ("error.txt","w");
    fprintf (pFile, "l2=%.4e,h1 = %.4e, surf = %.4e, dg = %.4e, sd = %.4e, total= %.4e,  ndof = %d, diam =%.4e",
Hermes::sqrt(err_l2_2),Hermes::sqrt(err_H1_2), Hermes::sqrt(err_surf_2),Hermes::sqrt(err_DG_2),Hermes::sqrt(err_sd_2), total,ndof,diam_min);
fclose (pFile);   
*/
// Output solution in VTK format.
	Linearizer lin;
	bool mode_3D = true;

MeshFunctionSharedPtr<double> filter(new AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_new, u_prev_time)));
//fview.show(filter);
lin.save_solution_vtk(u_new, "sln.vtk", "solution", mode_3D);
//lin.save_solution_vtk(u_prev_time, "init.vtk", "solution", mode_3D);
//lin.save_solution_vtk(filter, "error.vtk" , "error", false);  




}
