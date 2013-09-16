#include "definitions.h"


const double EPS = 1.;

const double penalty_parameter = 1.;

enum DG_TYPE {Baumann_Oden,	IP,	NIPG, CG, NONE};
DG_TYPE type = NIPG;



CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,bool all, bool DG) : WeakForm<double>(1)
{
 this->set_ext(sln_prev_time);

  if(all)
  {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
	add_vector_form(new RHS(0));
  } 
   if(DG) 	
		add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));
		
 	add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));   
	add_vector_form_surf(new CustomVectorFormSurface(0) );
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
		result += wt[i] * (EPS* (u->dx[i]*v->dx[i]+u->dy[i]*v->dy[i]));
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



//----------DG-Matrix-Form------------
template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
	Real flux_u = Real(0);
	double diam = e->diam;


  for (int i = 0; i < n; i++) 
  {



		Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
		flux_u = (u->fn_central == NULL ? (u->dx_neighbor[i]*e->nx[i]+u->dy_neighbor[i]* e->ny[i])  : (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i]) );
		Real jump_u =(u->fn_central == NULL ? -u->val_neighbor[i] :u->val[i]);
		Real mid_u = (u->fn_central == NULL ? u->val_neighbor[i] :u->val[i]);
		Real mid_v_dx =(v->fn_central == NULL ? (v->dx_neighbor[i]*e->nx[i]+v->dy_neighbor[i]* e->ny[i])  : (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]) ); 
		Real jump_v_dx =(v->fn_central == NULL ? -(v->dx_neighbor[i]*e->nx[i]+v->dy_neighbor[i]* e->ny[i])  : (v->dx[i]* e->nx[i]+v->dy[i]* e->ny[i]) ); 
		Real jump_u_dx= (u->fn_central == NULL ? -(u->dx_neighbor[i]*e->nx[i]+u->dy_neighbor[i]*e->ny[i])  : (u->dx[i]* e->nx[i]+u->dy[i]* e->ny[i]) );
		Real mid_v =(v->fn_central == NULL ? v->val_neighbor[i] : v->val[i]);
		Real mid_u_dx = flux_u;
		Real u_v_dx =0.;

Real jump_u_x = (u->fn_central == NULL ?-u->dx_neighbor[i] :u->dx[i]);
Real jump_u_y = (u->fn_central == NULL ?-u->dy_neighbor[i] :u->dy[i]);
Real jump_v_x = (v->fn_central == NULL ?-v->dx_neighbor[i] :v->dx[i]);
Real jump_v_y = (v->fn_central == NULL ?-v->dy_neighbor[i] :v->dy[i]);




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
  return result;
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



//--------Boundary-Matrix-Form
double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  	double result = 0.;
	double diam = e->diam;
		for (int i = 0; i < n; i++)
		{		
		
			if(type == Baumann_Oden)
			{	
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
			}else if(type == CG)
			{
				result += wt[i]*v->val[i]*u->val[i]/diam;
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

//----------Dirichlet Vector surface Form
double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
 Func<double>* exact = ext[0];	
			double diam = e->diam;
   for (int i = 0; i < n; i++)
	{ 
				

			if(type == Baumann_Oden)
			{	
					result += wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
			}else if(type == IP)
			{	
					result -= wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
					result += wt[i]*v->val[i]*exact->val[i]/diam*penalty_parameter;	
			}else if(type == NIPG)
			{
					result += wt[i]* EPS*(v->dx[i]*e->nx[i]+v->dy[i]* e->ny[i])*exact->val[i];
					result += wt[i]*v->val[i]*exact->val[i]/diam*penalty_parameter;		
			}else if(type == CG)
			{
					result += wt[i]* EPS*(exact->dx[i]*e->nx[i]+exact->dy[i]* e->ny[i])*v->val[i];
					result += wt[i]*v->val[i]*exact->val[i]/diam;		
			}

	}
  
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(10);
  //for (int i = 0; i < n; i++)
    //result += -wt[i] * v->val[i];
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
   for (int i = 0; i < n; i++)
		{ 	 double x= e->x[i]; double y = e->y[i];
			result += wt[i]*(2.*(2.-x*x-y*y))*v->val[i];
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

//-----------------error calculation-------------------------

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

 void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const {

			
		dx =(y*y-1.)*2.*x;
		dy =(x*x-1.)*2.*y;


};

 double CustomInitialCondition::value(double x, double y) const { 
       
  double result = 0.0;

result = (x*x-1.)*(y*y-1.);
	
return result;


};

 Ord CustomInitialCondition::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition::clone() const
    {
return new CustomInitialCondition(this->mesh);

    }
//------------------Boundary Condition------------




    double EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
{

		return exact_solution->value(x,y);

}




