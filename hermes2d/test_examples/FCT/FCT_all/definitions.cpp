#include "definitions.h"
const double EPS = 1e-6;

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

		//double diam = e->get_diameter();
		//double area =Hermes::sqrt(e->get_area());

		double abs_v =0;	
		double* x_coord = rm->get_phys_x(order);
		double* y_coord = rm->get_phys_y(order);
			for( int j = 0; j < np; ++j )
			{
				double v_x =2.; 
		double v_y =3.;	 
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}

//---------------Massematrix-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step,MeshFunctionSharedPtr<double> sln_prev_time) : WeakForm<double>(1) {
 this->set_ext(sln_prev_time);
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

 
//---------------KonvektionDiffusionReaction-----------
//----------------------wf-

  CustomWeakFormConvection::CustomWeakFormConvection(bool conv, bool diff,bool reac) : WeakForm<double>(1) {
    if(conv) add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
if(diff) add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0));
if(reac) add_matrix_form(new CustomMatrixFormVolReaction(0, 0));

  };

//------------Matrix Form  Konvektion------------
    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{		

Real v_x =Real(2.); 
Real v_y =Real(3.);

//result += wt[i]*u->val[i]*(v->dx[i]*v_x+ v->dy[i]*v_y); //mit surface
result -= wt[i]*v->val[i]*(u->dx[i]*v_x+ u->dy[i]*v_y); 

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



//------------Matrix Form //Diffusion------------
    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolDiffusion::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{	
result -= wt[i]*(v->dx[i]*u->dx[i]+ v->dy[i]*u->dy[i])*EPS; 

}
  return result;

    };

double CustomMatrixFormVolDiffusion::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomMatrixFormVolDiffusion::clone() const
{
  return new CustomMatrixFormVolDiffusion(*this);
}
//------------Matrix Form //REaction------------
    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolReaction::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{	
result -=wt[i]*u->val[i]*v->val[i];

}
  return result;

    };

double CustomMatrixFormVolReaction::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolReaction::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomMatrixFormVolReaction::clone() const
{
  return new CustomMatrixFormVolReaction(*this);
}



//------------surface---------
CustomWeakFormSurface::CustomWeakFormSurface(MeshFunctionSharedPtr<double> sln_prev_time): WeakForm<double>(1) {
 this->set_ext(sln_prev_time);
	// add_matrix_form_surf(new CustomMatrixFormSurface(0, 0)); 
		add_vector_form_surf(new CustomVectorFormSurface(0) );
  };
//matrix form surface

double CustomWeakFormSurface::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
	for (int i = 0; i < n; i++)
	{

		double v_x =2.; 
		double v_y =3.;
   double a_dot_n = (v_x*e->nx[i]+v_y*e->ny[i]);
		result -= wt[i] * static_cast<CustomWeakFormSurface*>(wf)->upwind_flux(u->val[i], 0., a_dot_n) * v->val[i];
	}
		return result;

}

Ord CustomWeakFormSurface::CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i]*u->val[i];
  return result;
}

MatrixFormSurf<double>* CustomWeakFormSurface::CustomMatrixFormSurface::clone() const
{
  return new CustomWeakFormSurface::CustomMatrixFormSurface(*this);
}
//----------Dirichlet Vector surface Form
double CustomWeakFormSurface::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
 	Func<double>* exact = ext[0];
   for (int i = 0; i < n; i++)
	{ 
		 //inlet
		double v_x =2.; 
		double v_y =3.;
			 double a_dot_n = (v_x*e->nx[i]+v_y*e->ny[i]);
				
	//result -= wt[i] * static_cast<CustomWeakFormSurface*>(wf)->upwind_flux(0, exact->val[i], a_dot_n) * v->val[i];
	result += wt[i]*(e->nx[i]*exact->dx[i]+ e->ny[i]*exact->dy[i])* v->val[i];
	}
  
  return result;
}

Ord CustomWeakFormSurface::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(10);
  return result;
}

VectorFormSurf<double>* CustomWeakFormSurface::CustomVectorFormSurface::clone() const
{
  return new CustomWeakFormSurface::CustomVectorFormSurface(*this);
}

double CustomWeakFormSurface::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakFormSurface::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}



//---------------wf rhs--------

CustomWeakFormRHS::CustomWeakFormRHS(MeshFunctionSharedPtr<double> sln_exact,MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,double time_step, double theta) : WeakForm<double>(1)
{
 this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(sln_exact,sln_prev_time) );
add_vector_form(new RHS(0, time_step, theta));

}

WeakForm<double>* CustomWeakFormRHS::clone() const
{
  return new CustomWeakFormRHS(*this);
}

//--------rhs Vector Form-----------------
double  CustomWeakFormRHS::RHS::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{  
	double result = 0; 
 Func<double>* exact = ext[0];
double time = static_cast<CustomWeakFormRHS*>(wf)->get_current_time(); 


   for (int i = 0; i < n; i++)
	{ 
		double v_x =2.; 
		double v_y =3.;
	
		//calculating f:
		double x = e->x[i]; double y = e->y[i];		
	double arg = 2./Hermes::sqrt(EPS)*(Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5));
	double u_t= 16.*Hermes::cos(PI*time)*PI*x*(1-x)*y*(1-y)*(0.5+Hermes::atan(arg)/PI);
		
		result +=wt[i] *(v->val[i]* (u_t+exact->val[i] + v_x*exact->dx[i] + v_y*exact->dy[i])+ EPS*(exact->dx[i]*v->dx[i]+exact->dy[i]*v->dy[i]));

	}
 return result;
}

    Ord  CustomWeakFormRHS::RHS::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
	return Ord(10);
}

    VectorFormVol<double>*  CustomWeakFormRHS::RHS::clone() const
{

	return new CustomWeakFormRHS::RHS(*this);
}


//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
	double t = time;
	double arg = 2./Hermes::sqrt(EPS)*(Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5));
	double c = 16.*Hermes::sin(PI*t);
	dx = c*y*(1.-y)*((1.-2.*x)*(0.5+Hermes::atan(arg)/PI)+(x-x*x)/PI*1./(1+arg*arg)*(-4*(x-0.5))/Hermes::sqrt(EPS));
	dy = c*x*(1.-x)*((1.-2.*y)*(0.5+Hermes::atan(arg)/PI)+(y-y*y)/PI*1./(1+arg*arg)*(-4*(y-0.5))/Hermes::sqrt(EPS));
};




 double CustomInitialCondition::value(double x, double y) const 
{
	double t = time;
	double result = 0.0;
	double arg = 2./Hermes::sqrt(EPS)*(Hermes::sqr(0.25)-Hermes::sqr(x-0.5)-Hermes::sqr(y-0.5));
	result= 16.*Hermes::sin(PI*t)*x*(1-x)*y*(1-y)*(0.5+Hermes::atan(arg)/PI);
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
