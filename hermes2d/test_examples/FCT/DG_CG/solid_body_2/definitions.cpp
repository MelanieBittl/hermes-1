#include "definitions.h"

const double EPS = 1e-3;
const double penalty_parameter = 1.;

enum DG_TYPE {Baumann_Oden,	IP,	NIPG, NONE};
DG_TYPE type =NIPG;


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

				double v_x = -y_coord[j];
				double v_y = x_coord[j];	 
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}

CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_exact,MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,double time_step, double theta,double theta_DG, bool all, bool DG, bool right_hand_side) : WeakForm<double>(1)
{
 this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(sln_exact,sln_prev_time) );

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
			 add_vector_form_DG(new CustomVectorFormInterface(0, theta_DG));
	}
		
   if(DG)				
		 add_matrix_form_DG(new CustomMatrixFormInterface(0, 0, theta_DG));	

}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}

//---------------volume integral-----------

    template<typename Real, typename Scalar>
    Scalar CustomWeakForm::CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{
	Real v_x =-e->y[i];
	Real v_y =e->x[i];
		result += wt[i] *( (u->val[i] *(v->val[i]/time_step - theta*(v->dx[i]*v_x+ v->dy[i]*v_y)))+ EPS*theta* (u->dx[i]*v->dx[i]+u->dy[i]*v->dy[i]) );
/*
double v_x =2.; double v_y =3.;
		result += wt[i] *( (u->val[i] *(v->val[i]/time_step + theta*v->val[i] - theta*(v->dx[i]*v_x+ v->dy[i]*v_y)))+ EPS*theta* (u->dx[i]*v->dx[i]+u->dy[i]*v->dy[i]) );*/
}
  return result;

    };

double CustomWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomWeakForm::CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {			
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}

//--------Boundary-Matrix-Form----
double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
	double diam = e->diam;
	for (int i = 0; i < n; i++)
	{
			double v_x =-e->y[i];	double v_y =e->x[i];
//double v_x =2.; double v_y =3.;

   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], 0., a_dot_n) * v->val[i];

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
template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
Real flux_u = Real(0);
	double diam = e->diam;


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
  return (result*theta);

}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
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
DiscontinuousFunc<double>* sln_prev_time = ext[1];		
  double result = double(0);
	double diam = e->diam;

  for (int i = 0; i < n; i++) 
  {
			double v_x =  (- e->y[i]);	double v_y = (e->x[i]) ;

//double v_x =2.; double v_y =3.; 

   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    double jump_v =  v->val[i];
	double	flux_u =(sln_prev_time->dx_neighbor[i]*e->nx[i]+sln_prev_time->dy_neighbor[i]* e->ny[i])  + (sln_prev_time->dx[i]* e->nx[i]+sln_prev_time->dy[i]* e->ny[i]);
	double jump_u =(sln_prev_time->val[i] -sln_prev_time->val_neighbor[i] );
	double mid_v_dx = (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]) ; 
double jump_v_dx = mid_v_dx;
  double mid_u = (sln_prev_time->val[i] + sln_prev_time->val_neighbor[i] );

double mid_u_dx = flux_u;
double jump_u_dx =  (sln_prev_time->dx[i]* e->nx[i]+sln_prev_time->dy[i]* e->ny[i])-(sln_prev_time->dx_neighbor[i]*e->nx[i]+sln_prev_time->dy_neighbor[i]* e->ny[i])  ;
double mid_v= v->val[i];

      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(sln_prev_time->val[i], sln_prev_time->val_neighbor[i], a_dot_n) * jump_v;	

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

//--------rhs Vector Form-----------------
double  CustomWeakForm::RHS::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{  
	double result = 0;
 	Func<double>* sln_prev_time = ext[1];
 Func<double>* exact = ext[0];
double t = static_cast<CustomWeakForm*>(wf)->get_current_time(); 


   for (int i = 0; i < n; i++)
		{ 

			double v_x =(- e->y[i]);			double v_y = (e->x[i]) ; 
			result += wt[i] *(sln_prev_time->val[i]* (v->val[i]/time_step +(1.-theta)*( v_x*v->dx[i] + v_y*v->dy[i]))- (1.-theta)*EPS*(sln_prev_time->dx[i]*v->dx[i]+sln_prev_time->dy[i]*v->dy[i]));

/*
double x = e->x[i]; double y = e->y[i];
double v_x =2.; double v_y =3.;

result += wt[i] *(sln_prev_time->val[i]* (v->val[i]/time_step- (1.-theta)*v->val[i] +(1.-theta)*( v_x*v->dx[i] + v_y*v->dy[i]))- (1.-theta)*EPS*(sln_prev_time->dx[i]*v->dx[i]+sln_prev_time->dy[i]*v->dy[i]));

//Berechnung von f
double a = 16.*Hermes::sin(M_PI*t)*x*(1.-x)*y*(1.-y);
double arg = Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5);
arg *= 2./Hermes::sqrt(EPS);
double b = (0.5+std::atan(arg)/M_PI);

double ax = 16.*Hermes::sin(M_PI*t)*y*(1.-y)*(1.-2.*x);
double axx = 16.*Hermes::sin(M_PI*t)*y*(1.-y)*(-2.);
double ay =	16.*Hermes::sin(M_PI*t)*x*(1.-x)*(1.-2.*y);
double ayy =	16.*Hermes::sin(M_PI*t)*x*(1.-x)*(-2.);

double bx = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-x+0.5);
double bxx = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-1.) - 1./(M_PI*Hermes::sqr(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-x+0.5)*2.*arg;
double by = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-y+0.5); 
double byy = 1./(M_PI*(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-1.) - 1./(M_PI*Hermes::sqr(1.+arg*arg))*2./Hermes::sqrt(EPS)*2.*(-y+0.5)*2.*arg;

double a_t = 16.*Hermes::cos(M_PI*t)*x*(1.-x)*y*(1.-y)*M_PI;
double laplace = axx*b+bxx*a+ax*bx*2. + ayy*b+byy*a+ay*by*2.;

result += wt[i]*(a_t*b+ 2.*exact->dx[i] + 3.*exact->dy[i] + exact->val[i]- laplace*EPS)*v->val[i] ;*/
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

CustomNormFormSurf_1::CustomNormFormSurf_1(int i, int j) : NormFormSurf<double>(i, j)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormDG_1::CustomNormFormDG_1(int i, int j) : NormFormDG<double>(i, j)
{
}

CustomNormFormSurf_2::CustomNormFormSurf_2(int i, int j) : NormFormSurf<double>(i, j)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormDG_2::CustomNormFormDG_2(int i, int j) : NormFormDG<double>(i, j)
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
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem); 

   double result = double(0);
  for (int i = 0; i < n; i++)
	{			
			double v_x =(- e->y[i]);
			double v_y = (e->x[i]) ; 
    result += wt[i] * (v_x*u->dx[i] + v_y*u->dy[i]) * (v_x*v->dx[i] + v_y*v->dy[i]);
	}
  return (result*diam/(abs_v));


}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
			double v_x =(- e->y[i]);
			double v_y = (e->x[i]) ; 
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
			
			double v_x =(- e->y[i]);
			double v_y = (e->x[i]) ; 
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * a_dot_n;
		}
  return result;
}

/////////////
double CustomNormFormSurf_1::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
				
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
    result += wt[i] * u->val[i] * v->val[i] * penalty_parameter;
	}
  return result;
}

double CustomNormFormDG_1::value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
{
 
  double result = double(0);
  for (int i = 0; i < n; i++)
	{

    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * penalty_parameter;
		}
  return result;
}
double CustomNormFormSurf_2::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
	double mid_v_dx = (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]); 
	double mid_u_dx = (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i]); 
    result += wt[i] * mid_v_dx* mid_u_dx * 1./penalty_parameter;
	}
  return result;
}

double CustomNormFormDG_2::value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
{
 
  double result = double(0);
  for (int i = 0; i < n; i++)
	{
	double mid_v_dx = (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i])+(v->dx_neighbor[i]*e->nx[i]+v->dy_neighbor[i]* e->ny[i]); 
	double mid_u_dx = (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i])+(u->dx_neighbor[i]*e->nx[i]+u->dy_neighbor[i]* e->ny[i]); 
    result += wt[i] * mid_v_dx/2.*mid_u_dx/2. * 1./penalty_parameter;
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

