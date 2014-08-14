#include "definitions.h"
double calc_abs_v(Element* e, double time)
{
	double final_time = 1.5;
Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;
			if(e->is_triangle()) mode = HERMES_MODE_TRIANGLE;
	// order of integral
	const int order = 4;
	// refmap for computing Jacobian
	RefMap* rm = new RefMap;
	rm->set_quad_2d(&g_quad_2d_std);

			// get the quadrature points
			int np = g_quad_2d_std.get_num_points(order,mode);
			double3 *pt = g_quad_2d_std.get_points(order, mode);
			// get the constant Jacobian
		  rm->set_active_element(e);

		double diam = e->get_diameter();
		double area =Hermes::sqrt(e->get_area());

		double abs_v =0;	
		double* x_coord = rm->get_phys_x(order);
		double* y_coord = rm->get_phys_y(order);
			for( int j = 0; j < np; ++j )
			{	double x = x_coord[j]; double y = y_coord[j];
				double g = Hermes::cos(PI*time/final_time);
				double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
				double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}

CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,double final_time,double time_step, double theta, bool all, bool DG, bool SD,bool right_hand_side) : WeakForm<double>(1), final_time(final_time)
{
 this->set_ext(sln_prev_time);

  if(all) 
	{ 
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0, time_step, theta));
	 add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));   

	}
  
	if(right_hand_side)
	{
		add_vector_form(new RHS(0, time_step, theta));
		add_vector_form_surf(new CustomVectorFormSurface(0) );
			if(DG)
		 add_vector_form_DG(new CustomVectorFormInterface(0, theta));
	}
		
   if(DG)			
		 add_matrix_form_DG(new CustomMatrixFormInterface(0, 0,theta));


		
		
	 if(SD)	add_matrix_form(new Streamline(0,0,mesh));
   

}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}

//---------------mass+convection-----------


double CustomWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
          double result = 0.; 
	  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
  for (int i = 0; i < n; i++)
{
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;

		result += wt[i] * (u->val[i] *(v->val[i]/time_step - theta*(v->dx[i]*v_x+ v->dy[i]*v_y)));

}
  return result;
    };

   Ord CustomWeakForm::CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {			
     return Ord(10);
    };

   MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}

//---------Streamline Diffusion---------------

double CustomWeakForm::Streamline::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
   double result = 0.;
double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
		double diam = e->diam;
	
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem,time);
for (int i = 0; i < n; i++)
{
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
double tau = diam/(2.*abs_v);
  result += wt[i] *(tau*(u->dx[i] * v_x + u->dy[i] * v_y ) *(v->dx[i] * v_x + v->dy[i] * v_y ));

}


  return result;
}

Ord CustomWeakForm::Streamline::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
 return Ord(10);
}

MatrixFormVol<double>* CustomWeakForm::Streamline::clone() const
{
  return new CustomWeakForm::Streamline(*this);
}




//--------Boundary-Matrix-Form-----outlet
double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
	for (int i = 0; i < n; i++)
	{
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] *static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], 0., a_dot_n) * v->val[i];
		
	}
		return result;

}

Ord CustomWeakForm::CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i]*u->val[i];
  return result;
}

MatrixFormSurf<double>* CustomWeakForm::CustomMatrixFormSurface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormSurface(*this);
}
//----------DG-Matrix-Form------------
double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
  for (int i = 0; i < n; i++) 
  {
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0., u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i],0., a_dot_n) * jump_v;
      
  }
  return result*theta;
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return Ord(10);
}

MatrixFormDG<double>* CustomWeakForm::CustomMatrixFormInterface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormInterface(*this);
}



//------------DG Vector Form--------------
     double CustomWeakForm::CustomVectorFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
        Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
DiscontinuousFunc<double>* exact = ext[0];		
  double result = double(0);
	double diam = e->diam;
	  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
  for (int i = 0; i < n; i++) 
  {
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    double jump_v =  v->val[i];
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(exact->val[i], exact->val_neighbor[i], a_dot_n) * jump_v;
      
  }
  return (-result*(1.-theta));


}

    Ord CustomWeakForm::CustomVectorFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const
{  
Ord result = Ord(10); 
  return result;
}

      VectorFormDG<double>* CustomWeakForm::CustomVectorFormInterface::clone() const
{
  return new CustomWeakForm::CustomVectorFormInterface(*this);
}

//----------Dirichlet Vector surface Form
double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
 Func<double>* exact = ext[0];	
 	  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
   for (int i = 0; i < n; i++)
	{ 
		 //inlet
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
			
			double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
				
			result -= wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0, exact->val[i], a_dot_n) * v->val[i];

	}
  
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(10);
  return result;
}

VectorFormSurf<double>* CustomWeakForm::CustomVectorFormSurface::clone() const
{
  return new CustomWeakForm::CustomVectorFormSurface(*this);
}

//--------rhs Vector Form-----------------
double  CustomWeakForm::RHS::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{  
	double result = 0;
 	Func<double>* exact = ext[0];
	 	  double time = static_cast<CustomWeakForm*>(wf)->get_current_time(); 
	  double final_time = static_cast<CustomWeakForm*>(wf)->get_final_time(); 
   for (int i = 0; i < n; i++)
		{ 
		double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
			result += wt[i] *exact->val[i]* (v->val[i]/time_step +(1.-theta)*( v_x*v->dx[i] + v_y*v->dy[i])) ;
		}
 return result;
}

    Ord  CustomWeakForm::RHS::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
	return Ord(10);
}

    VectorFormVol<double>*  CustomWeakForm::RHS::clone() const
{

	return new CustomWeakForm::RHS(*this);
}

///----------------helper functions--------------
double CustomWeakForm::calculate_a_dot_v(double x, double y, double vx, double vy) const
{
 
 return  x*vx + y*vy;
}

Ord CustomWeakForm::calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const
{
  return Ord(10);
}

double CustomWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}


//---------------Massematrix-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix() : WeakForm<double>(1) {
		add_matrix_form(new DefaultMatrixFormVol<double>(0,0)); 	

  }

	    WeakForm<double>* CustomWeakFormMassmatrix::clone() const
    {
    CustomWeakFormMassmatrix* wf;
    wf = new CustomWeakFormMassmatrix;

    return wf;
    }

//--------------error_calculation----------------------


CustomNormFormVol::CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
{
  this->set_area(HERMES_ANY);
}
StreamlineDiffusionNorm::StreamlineDiffusionNorm(int i, int j,MeshSharedPtr mesh, double final_time, double time) : NormFormVol<double>(i, j), mesh(mesh), time(time), final_time(final_time)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormSurf::CustomNormFormSurf(int i, int j, double final_time, double time) : NormFormSurf<double>(i, j), time(time), final_time(final_time)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormDG::CustomNormFormDG(int i, int j, double final_time, double time) : NormFormDG<double>(i, j), time(time), final_time(final_time)
{
}

double CustomNormFormVol::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i];
  return result;
}



double StreamlineDiffusionNorm::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
		double diam = e->diam;
		//double area = Hermes::sqrt(e->area);
		Element* elem = mesh->get_element(e->id);		
		double abs_v =calc_abs_v(elem,time); 

   double result = double(0);
  for (int i = 0; i < n; i++)
	{				double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
    result += wt[i] * (v_x*u->dx[i] + v_y*u->dy[i]) * (v_x*v->dx[i] + v_y*v->dy[i]);
	}
  return (result*diam/(abs_v));
		

}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
			double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * u->val[i] * v->val[i] * a_dot_n;
	}
  return result;
}

double CustomNormFormDG::value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
{
 
  double result = double(0);
  for (int i = 0; i < n; i++)
	{
			double x = e->x[i]; double y = e->y[i];
		double g = Hermes::cos(PI*time/final_time);
		double v_x = Hermes::sin(PI*x)*Hermes::sin(PI*x)*Hermes::sin(2.*PI*y)*g;
		double v_y = -Hermes::sin(PI*y)*Hermes::sin(PI*y)*Hermes::sin(2.*PI*x)*g;
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * a_dot_n;
		}
  return result;
}




//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
		dx=0.; dy=0.;
}

 double CustomInitialCondition::value(double x, double y) const 
{
       
  double result = 0.0;
	double radius_2;
        //hump
	double x_0 =1.;
	double y_0= 1.;	
	radius_2 =  std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0);
	if( radius_2<= 0.8) 	
		 result = 1.;
       return result;

}

 Ord CustomInitialCondition::ord(double x, double y) const 
 {
      return Ord(10);
	}

 MeshFunction<double>* CustomInitialCondition::clone() const
	{
		return new CustomInitialCondition(this->mesh);

	}



