#include "definitions.h"

CustomWeakForm::CustomWeakForm(double time_step,double theta, Solution<double>* sln_prev_time, std::string inlet, Mesh* mesh, bool all, bool DG) : WeakForm<double>(1), mesh(mesh)
{
 this->set_ext(sln_prev_time);

  if(all)
  {
    add_matrix_form(new CustomMatrixFormVol(0, 0, time_step,theta));
    }
   add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));    
   if(DG) add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));
   
	add_vector_form_surf(new CustomVectorFormSurface(0,inlet) );
}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
{
			Real v_x = e->y[i]; 
			Real v_y = Real(1.)-e->x[i];
    result -= wt[i] *(u->val[i] *(v->dx[i] * v_x + v->dy[i] * v_y ));
	//result += wt[i] *(v->val[i] *(u->dx[i] * v_x + u->dy[i] * v_y ));

}


  return result;
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVol::clone() const
{
  return new CustomWeakForm::CustomMatrixFormVol(*this);
}




double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
 
	for (int i = 0; i < n; i++)
	{

	//Ausstroemrand
if((e->y[i]==1)||(e->x[i]==1))
	{	
	double v_x = e->y[i];
 	double v_y = 1.-e->x[i]; 
   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] *static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], 0., a_dot_n) * v->val[i];	
	//result += wt[i] *u->val[i] *a_dot_n* v->val[i];	
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

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                        Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
 	 Real v_x = (e->y[i]);
 	 Real v_y = Real(1.)-e->x[i]; 
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
   result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v; 
      
  }
  return result;
}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                             Geom<Ord> *e, Func<Ord> **ext) const
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
   for (int i = 0; i < n; i++){ 
 //Dirichlet-Rand!
if((e->x[i]<1)&&(e->y[i]==0)) 
		{
			double v_x =	e->y[i]; 
			double v_y = 1.-e->x[i];
			double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
			result -= wt[i] * exact->val[i] * v->val[i] * a_dot_n;
		}
	}
  
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += -wt[i] * v->val[i];
  return result;
}

VectorFormSurf<double>* CustomWeakForm::CustomVectorFormSurface::clone() const
{
  return new CustomWeakForm::CustomVectorFormSurface(*this);
}

template<typename Real>
Real CustomWeakForm::CustomVectorFormSurface::F(Real x, Real y) const
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormSurface::g(std::string ess_bdy_marker) const
{
  if(ess_bdy_marker == inlet) return 1; else return 0;
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



//---------------Konvektion-----------

  CustomWeakFormConvection::CustomWeakFormConvection() : WeakForm<double>(1) 
{
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));

  };



    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
  	result -= wt[i] * (u->val[i] *(v->dx[i] * (e->y[i]) + v->dy[i] * (1.-e->x[i]) ));
	//result -= wt[i] * (u->val[i] *(v->dx[i]*0.5 + v->dy[i]  ));
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







//------------------------------Boundary Condition------------------------
CustomDirichletCondition::CustomDirichletCondition(Hermes::vector<std::string> markers)
  : EssentialBoundaryCondition<double>(markers)
{
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{
  return EssentialBoundaryCondition<double>::BC_FUNCTION;
}

double CustomDirichletCondition::value(double x, double y, double n_x, double n_y,
                                       double t_x, double t_y) const
{
	double result = 0;
		double radius = Hermes::sqrt((x-1)*(x-1)+y*y);
		if((radius>= 0.5)&&(radius<=0.8))
		{		
			double arg = PI*(radius-0.65)/0.15;
			result = 0.25*(1+Hermes::cos(arg));
		}	
		
	return result;
}

//------------------- Initial condition ----------------

 void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const {
      
			double radius = Hermes::sqrt((x-1)*(x-1)+y*y);
	/*if((radius>= 0.5)&&(radius<=0.8)){		
		double arg = PI*(radius-0.65)/0.15;
		dx = -0.25*Hermes::sin(arg)*(PI/0.15)*x/radius;
		dy	= -0.25*Hermes::sin(arg)*(PI/0.15)*y/radius;
	}*/
	if((radius>= 0.1)&&(radius<=0.9))
	{		
		double arg = PI*(radius-0.5)/0.8;
		dx = -Hermes::sin(arg)*(PI/0.3)*x/radius;
		dy	= -Hermes::sin(arg)*(PI/0.3)*y/radius;
		double result= Hermes::cos(arg);
		int k = 4;
		dx *=k*std::pow(result, k-1);
		dy *=k*std::pow(result, k-1);
	}

	
// fuer v = (0.5,1)
/*	if((x-0.5*y>0)&&(x-0.5*y<0.5))
{
	double arg = PI*(x-0.5*y-0.25)/0.5;
	dx= -std::sin(arg)*PI/0.5;
	dy=	std::sin(arg)*PI/0.5;
	int k =2;
	dx*= k*std::pow(std::cos(arg),k-1);
	dy *=k*std::pow(std::cos(arg),k-1);
}
// fuer v = (1,1)
/*if((x-y>0)&&(x-y<0.5))
{
	double arg = PI*(x-y-0.25)/0.5;
	dx= -std::sin(arg)*PI/0.5;
	dy=	std::sin(arg)*PI/0.5;
	int k =2;
	dx*= k*std::pow(std::cos(arg),k-1);
	dy *=k*std::pow(std::cos(arg),k-1);
}*/
		else 
	{
		dx =0;
		dy =0;
	}	

};

 double CustomInitialCondition::value(double x, double y) const {
       
 		 double result = 0.0;
		double radius = Hermes::sqrt((x-1)*(x-1)+y*y);
		/*if((radius>= 0.5)&&(radius<=0.8))
		{		
			double arg = PI*(radius-0.65)/0.15;
			result = 0.25*(1+Hermes::cos(arg));
		}	
		
*/
	if((radius> 0.1)&&(radius<0.9))
		{
		
		double arg = PI*(radius-0.5)/0.8;
		result= Hermes::cos(arg);
		}
		int k = 4;
		return std::pow(result,k);



// fuer v = (0.5,1)
/*	if((x-0.5*y>0)&&(x-0.5*y<0.5))
	{
		double arg = PI*(x-0.5*y-0.25)/0.5;
		result= std::cos(arg);
		int k = 2;
		return std::pow(result,k);		
	}

// fuer v = (1,1)
/*if((x-y>0)&&(x-y<0.5))
	{
		double arg = PI*(x-y-0.25)/0.5;
		result= std::cos(arg);
		int k = 2;
		return std::pow(result,k);		
	}*/
	
	return result;


};

 Ord CustomInitialCondition::ord(Ord x, Ord y) const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition::clone() const
    {
return new CustomInitialCondition(this->mesh);

    }


 void CustomV_x::derivatives(double x, double y, double& dx, double& dy) const 
 {   
	dx=0.; dy=1.0;		
};

 double CustomV_x::value(double x, double y) const 
 {
       return (y);
};

 Ord CustomV_x::ord(Ord x, Ord y) const {
      return Ord(2);
};

 MeshFunction<double>* CustomV_x::clone() const
    {
return new CustomV_x(this->mesh);

    }
 void CustomV_y::derivatives(double x, double y, double& dx, double& dy) const 
 {   
	dx=-1.; dy=0.0;		

};

 double CustomV_y::value(double x, double y) const 
 {       
		return (-x);
};

 Ord CustomV_y::ord(Ord x, Ord y) const {
      return Ord(2);
};

 MeshFunction<double>* CustomV_y::clone() const
    {
return new CustomV_y(this->mesh);

    }
