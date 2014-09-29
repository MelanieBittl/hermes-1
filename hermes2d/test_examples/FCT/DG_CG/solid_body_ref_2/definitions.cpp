#include "definitions.h"

const double EPS = 1e-3;
const double penalty_parameter = 1.;

enum DG_TYPE {Baumann_Oden,	IP,	NIPG, NONE};
DG_TYPE type = Baumann_Oden;

double calc_abs_v(Element* e)
{
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
			{
				double v_x = y_coord[j]; 
				double v_y = 1.-x_coord[j];	
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}


CustomWeakForm::CustomWeakForm(bool DG) : WeakForm<double>(1), mesh(mesh)
{

   add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));    
   if(DG) add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));

}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}



template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormSurface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                      Geom<Real> *e, Func<Scalar> **ext) const
{

 	Scalar result = Scalar(0);
	Scalar diam = e->diam;
	for (int i = 0; i < n; i++)
	{
			Real  v_x =-e->y[i];	Real v_y =e->x[i];
//Real v_x =Real(2.); Real v_y =Real(3.);

   Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i],  Scalar(0), a_dot_n) * v->val[i];

			if(type == Baumann_Oden)
			{					result += wt[i]*v->val[i]*u->val[i]/diam*penalty_parameter;
				result += wt[i]*EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i]) *u->val[i];
				result -= wt[i]*EPS*(u->dx[i]*e->nx[i]+u->dy[i]* e->ny[i]) *v->val[i];
			}else if(type == IP)
			{	
				result -= wt[i]*EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i]) *u->val[i];
				result -= wt[i]*EPS*(u->dx[i]*e->nx[i]+u->dy[i]* e->ny[i]) *v->val[i];
				result += wt[i]*v->val[i]*u->val[i]/diam*penalty_parameter;
			}else if(type == NIPG)
			{
				result += wt[i]*EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i]) *u->val[i];
				result -= wt[i]*EPS*(u->dx[i]*e->nx[i]+u->dy[i]* e->ny[i]) *v->val[i];
				result += wt[i]*v->val[i]*u->val[i]/diam*penalty_parameter;
			}
	}
	return result;

}

double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* CustomWeakForm::CustomMatrixFormSurface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormSurface(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
Real flux_u = Real(0);
	Scalar diam = e->diam;


  for (int i = 0; i < n; i++) 
  {
	Real v_x =-e->y[i];	Real v_y =e->x[i];

//double v_x =2.; double v_y =3.;

    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);

    Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
		flux_u = (u->fn_central == NULL ? (u->dx_neighbor[i]*e->nx[i]+u->dy_neighbor[i]* e->ny[i])  : (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i]) );
	Real jump_u =(u->fn_central == NULL ? -u->val_neighbor[i] :u->val[i]);
	Real mid_u = (u->fn_central == NULL ? u->val_neighbor[i] :u->val[i]);
	Real mid_v_dx =(v->fn_central == NULL ? (v->dx_neighbor[i]*e->nx[i]+v->dy_neighbor[i]* e->ny[i])  : (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]) ); 
	Real jump_v_dx =(v->fn_central == NULL ? -(v->dx_neighbor[i]*e->nx[i]+v->dy_neighbor[i]* e->ny[i])  : (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]) ); 
Real jump_u_dx= (u->fn_central == NULL ? -(u->dx_neighbor[i]*e->nx[i]+u->dy_neighbor[i]* e->ny[i])  : (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i]) );
Real mid_v =(v->fn_central == NULL ? v->val_neighbor[i] : v->val[i]);

Real mid_u_dx = flux_u;

		if(type == Baumann_Oden)
		{
			result -= wt[i]*EPS*flux_u*jump_v/2.;
			result += wt[i]*EPS*jump_u*mid_v_dx/2.;
		}else if(type == IP)
		{
			result -= wt[i]*EPS*flux_u*jump_v/2.;
			result -= wt[i]*EPS*jump_u*mid_v_dx/2.;
			result += wt[i]*jump_u/diam*jump_v*penalty_parameter;
		}else if(type == NIPG)
		{
			result -= wt[i]*EPS*flux_u*jump_v/2.;
			result += wt[i]*EPS*jump_u*mid_v_dx/2.;
			result += wt[i]*jump_u/diam*jump_v*penalty_parameter;
		}



    if(u->fn_central == NULL)
		{
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(Real(0), u->val_neighbor[i], a_dot_n)* jump_v;
    }else{
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Real(0), a_dot_n) * jump_v;
		}
      
  }
  return (result);
}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormDG<double>* CustomWeakForm::CustomMatrixFormInterface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormInterface(*this);
}

double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
 Func<double>* exact = ext[0];
	double diam = e->diam;
   for (int i = 0; i < n; i++)
	{ 		 
			double v_x = -e->y[i];			double v_y = e->x[i]; 
//double v_x =2.; double v_y =3.;

			double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
			double grad_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(exact->dx[i], exact->dy[i], e->nx[i], e->ny[i]);


				result -= wt[i]*(static_cast<CustomWeakForm*>(wf)->upwind_flux(0., exact->val[i], a_dot_n))*v->val[i];

			if(type == Baumann_Oden)
			{	result += wt[i]*v->val[i]*exact->val[i]/diam*penalty_parameter;
					result += wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
			}else if(type == IP)
			{	
					result -= wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
					result += wt[i]*v->val[i]*exact->val[i]/diam*penalty_parameter;	
			}else if(type == NIPG)
			{
					result += wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
					result += wt[i]*v->val[i]*exact->val[i]/diam*penalty_parameter;		
			}

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





//---------------mass-matrix/tau-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step) : WeakForm<double>(1) 
 {
		CustomMatrixFormVolMassmatrix* mass_form= new CustomMatrixFormVolMassmatrix(0, 0, time_step);	
		add_matrix_form(mass_form);	 

  }

    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolMassmatrix::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
		     Scalar result = Scalar(0); 
	  for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->val[i])/time_step;
	  return result;

    };

   double CustomMatrixFormVolMassmatrix::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

    Ord CustomMatrixFormVolMassmatrix::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

MatrixFormVol<double>* CustomMatrixFormVolMassmatrix::clone() const
{
  return new CustomMatrixFormVolMassmatrix(*this);
}



//---------------Convection-----------

 
  CustomWeakFormConvection::CustomWeakFormConvection() : WeakForm<double>(1) 
  {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  };



    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{
Real  v_x =-e->y[i];	Real v_y =e->x[i];
result += wt[i]*u->val[i]*(v->dx[i]*v_x + v->dy[i]*v_y); 
result -= wt[i]*(v->dx[i]*u->dx[i]+ v->dy[i]*u->dy[i])*EPS; 
}
  return result;

    };

double CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}




//--------------error_calculation----------------------


CustomNormFormVol::CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
{
  this->set_area(HERMES_ANY);
}
StreamlineDiffusionNorm::StreamlineDiffusionNorm(int i, int j,MeshSharedPtr mesh) : NormFormVol<double>(i, j), mesh(mesh)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormSurf::CustomNormFormSurf(int i, int j) : NormFormSurf<double>(i, j)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormDG::CustomNormFormDG(int i, int j) : NormFormDG<double>(i, j)
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
		//double v_x = 0.5; 
	//	double v_y = 1.;
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem); //Hermes::sqrt(v_x*v_x+v_y*v_y);

   double result = double(0);
  for (int i = 0; i < n; i++)
	{			double v_x = -e->y[i];
 			double v_y = e->x[i]; 
    result += wt[i] * (v_x*u->dx[i] + v_y*u->dy[i]) * (v_x*v->dx[i] + v_y*v->dy[i]);
	}
  return (result*diam/(abs_v));
		//return result;

}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
	double v_x = -e->y[i];
 			double v_y = e->x[i]; 
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
			double v_x = -e->y[i];
 			double v_y = e->x[i];  
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * a_dot_n;
		}
  return result;
}


//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
	double t = time;
	double x_real = x_0*Hermes::cos(t) - y_0*Hermes::sin(t);
	double y_real = y_0*Hermes::cos(t) - x_0*Hermes::sin(t);
	double radius_2 = Hermes::sqr(x-x_real) + Hermes::sqr(y-y_real);

dx = -2*PI/Hermes::sqr(4.*PI*EPS*t)*std::exp(-radius_2/(4.*EPS*t))*(x-x_real);
dy = -2*PI/Hermes::sqr(4.*PI*EPS*t)*std::exp(-radius_2/(4.*EPS*t))*(y-y_real);

/*
double a = 16.*Hermes::sin(M_PI*t)*x*(1.-x)*y*(1.-y);
double arg = Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5);
arg *= 2./Hermes::sqrt(EPS);
double b = (0.5+std::atan(arg)/M_PI);

double ax = 16.*Hermes::sin(M_PI*t)*y*(1.-y)*(1.-2.*x);
double ay =	16.*Hermes::sin(M_PI*t)*x*(1.-x)*(1.-2.*y);

double bx = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-x+0.5);
double by = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-y+0.5); 
dx = ax*b + a*bx;
dy = ay*b + a*by;*/
};




 double CustomInitialCondition::value(double x, double y) const 
{
       double t = time;
  double result = 0.0;
	double x_real = x_0*Hermes::cos(t) - y_0*Hermes::sin(t);
	double y_real = y_0*Hermes::cos(t) - x_0*Hermes::sin(t);
	double radius_2 = Hermes::sqr(x-x_real) + Hermes::sqr(y-y_real);
		result= 1./(4.*PI*EPS*t)*std::exp(-radius_2/(4.*EPS*t));


/*
double a = 16.*Hermes::sin(M_PI*t)*x*(1.-x)*y*(1.-y);
double arg = Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5);
arg *= 2./Hermes::sqrt(EPS);
double b = (0.5+std::atan(arg)/M_PI);
result = a*b;*/

       return result;


};

 Ord CustomInitialCondition::ord(double x, double y)  const 
 {
      return Ord(10);
	};
 MeshFunction<double>* CustomInitialCondition::clone() const
	{
		return new CustomInitialCondition(this->mesh, this->time);

	}


