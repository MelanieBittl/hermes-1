#include "definitions.h"

CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_prev_time,bool all, bool DG) : WeakForm<double>(1)
{
 this->set_ext(sln_prev_time);

  if(all)
  {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
		//add_vector_form(new RHS(0));
    }
   add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));    
   if(DG) add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));
   
	add_vector_form_surf(new CustomVectorFormSurface(0) );
}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}



double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
double v_x =1.;
double v_y = 1.;
 
	for (int i = 0; i < n; i++)
	{

	//outlet
if((e->y[i]==1)||(e->x[i]==1))
	{	
	//double v_x = e->y[i];
 	//double v_y = 1.-e->x[i]; 
   double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] *static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], 0., a_dot_n) * v->val[i];
		}
	}
		return result;

}

Ord CustomWeakForm::CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(10);
  //for (int i = 0; i < n; i++)
    //result += wt[i] * v->val[i]*u->val[i];
  return result;
}

MatrixFormSurf<double>* CustomWeakForm::CustomMatrixFormSurface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormSurface(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
	Real v_x = Real(1);
	Real v_y = Real(1);
  for (int i = 0; i < n; i++) 
  {
 	 //Real v_x = (e->y[i]);
 	// Real v_y = Real(1.)-e->x[i]; 
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   // Real jump_v = v->val[i] - v->val_neighbor[i];
   //result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], u->val_neighbor[i], a_dot_n) * jump_v; 
    Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(Real(0), u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Real(0), a_dot_n) * jump_v;
      
  }
  return result;
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
	double v_x =1.;
	double v_y = 1.;		
   for (int i = 0; i < n; i++)
	{ 
		 //Dirichlet-Bdry!
		if((e->y[i]==0)||(e->x[i]==0))
				{
					//double v_x = e->y[i]; 
					//double v_y = 1.-e->x[i];
					double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
					result -= wt[i] * exact->val[i] * v->val[i] * a_dot_n;
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


double  CustomWeakForm::RHS::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{  
	double result = 0;
 	Func<double>* exact = ext[0];
   for (int i = 0; i < n; i++)
		{ 
			double x = e->x[i]; double y = e->y[i];
	/*if((x-0.5*y>0)&&(x-0.5*y<0.5))
		{
			double arg = PI*(x-0.5*y-0.25)*2.;
			double erg = (Hermes::cos(arg)-2*PI*Hermes::sin(arg))*Hermes::cos(arg);
			result += wt[i] * erg * v->val[i];
		}*/
			result += wt[i] * exact->val[i] * v->val[i];
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



//---------------convection-----------

  CustomWeakFormConvection::CustomWeakFormConvection() : WeakForm<double>(1) 
{
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));

  };



    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
		//result += wt[i] * (v->val[i] *((u->dx[i]+ u->dy[i])));
		result -= wt[i] * (u->val[i] *((v->dx[i]+ v->dy[i])));
		//result += wt[i] * (u->val[i] *(v->val[i]-(v->dx[i]+ v->dy[i])));
  //result -= wt[i] * (u->val[i] *(v->dx[i] * (e->y[i]) + v->dy[i] * (1.-e->x[i]) ));
  return result;

    };

double CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
			return Ord(10);
      //return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}


//--------------------REsidual--------------------------

  template<typename Real, typename Scalar>
  Scalar Residual_Mat::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const
{
     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
		result += wt[i] * (v->val[i] *((u->dx[i]+ u->dy[i])));
  return result;


}

 double Residual_Mat::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
 return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord Residual_Mat::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, Func<Ord>  **ext) const
{
return Ord(10);
}

    MatrixFormVol<double>* Residual_Mat::clone() const
{
	return new Residual_Mat(*this);
}



	double Residual::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{

	Func<double>* sln = ext[0];	Func<double>* exact = ext[1];
  double result = double(0);
  for (int i = 0; i < n; i++)
	{		//if((e->y[i]==0)||(e->x[i]==0))
			//result += wt[i] * (v->val[i] *(exact->dx[i]+exact->dy[i]));
	//else 
    	//result += wt[i] * (v->val[i] *(sln->dx[i]+sln->dy[i]));

	/*if((e->y[i]==0)||(e->x[i]==0))
			result -= wt[i] * (exact->val[i] *(v->dx[i]+v->dy[i]));
	else */
    	result -= wt[i] * (sln->val[i] *(v->dx[i]+v->dy[i]));
		}
  return result;

}

  Ord Residual::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{

return Ord(10);
}

 VectorFormVol<double>* Residual::clone() const
	{	return new Residual(*this);
		}


	double Residual_surf::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
	Func<double>* sln = ext[0];
	Func<double>* exact = ext[1];
double v_x =1.;
double v_y = 1.;		
   for (int i = 0; i < n; i++){ 
			double a_dot_n = v_x*e->nx[i]+ v_y*e->ny[i];
	//if((e->y[i]==0)||(e->x[i]==0))
	//	result += wt[i] * exact->val[i] * v->val[i] * a_dot_n;
	//else 
			result += wt[i] * sln->val[i] * v->val[i]* a_dot_n;
		
	}
  
  return result;

}

  Ord Residual_surf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{

return Ord(10);
}

 VectorFormSurf<double>* Residual_surf::clone() const
	{	return new Residual_surf(*this);
		}


  Wf_residual::Wf_residual(MeshFunctionSharedPtr<double> sln_1,MeshFunctionSharedPtr<double> sln_2 ): WeakForm<double>(1)
{ 
		this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(sln_1,sln_2));
		add_vector_form(new Residual(0));
		add_vector_form_surf(new Residual_surf(0));
		add_matrix_form(new Residual_Mat(0,0));
}

WeakForm<double>*   Wf_residual::clone() const
{
  return new   Wf_residual(*this);
}


//--------------error_calculation----------------------


CustomNormFormVol::CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
{
  this->set_area(HERMES_ANY);
}
StreamlineDiffusionNorm::StreamlineDiffusionNorm(int i, int j) : NormFormVol<double>(i, j)
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
		double area = Hermes::sqrt(e->area);
		double abs_v = Hermes::sqrt(2);
   double result = double(0);
  for (int i = 0; i < n; i++)
	{	
    result += wt[i] *(u->dx[i]+u->dy[i]) * (v->dx[i]+v->dy[i]) ;
	}
  return (result*diam/(area*abs_v));
		//return result;

}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
			double v_x = 1.; 
			double v_y = 1.;
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * u->val[i] * v->val[i]*a_dot_n;
	}
  return result;
}

double CustomNormFormDG::value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
{
 
  double result = double(0);
  for (int i = 0; i < n; i++)
	{
		double v_x = 1.; 
		double v_y = 1.;
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i])*a_dot_n;
		}
  return result;
}

//------------------- Initial condition ----------------

 void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const {
      
/*	double radius = Hermes::sqrt((x-1)*(x-1)+y*y);
	if((radius>= 0.1)&&(radius<=0.9))
	{		
		double arg = PI*(radius-0.5)/0.8;
		dx = -Hermes::sin(arg)*(PI/0.3)*x/radius;
		dy	= -Hermes::sin(arg)*(PI/0.3)*y/radius;
		double result= Hermes::cos(arg);
		int k = 4;
		dx *=k*std::pow(result, k-1);
		dy *=k*std::pow(result, k-1);
	}*/
	
	if((x-y>0)&&(x-y<0.5))
	{
		double arg = PI*(x-y-0.25)*2;
		dx= -std::sin(arg)*PI*2;
		dy=	std::sin(arg)*PI*2;
		int k =2;
		dx *=k*std::pow(std::cos(arg),k-1);
		dy *=k*std::pow(std::cos(arg),k-1);
	}
/*
if((x-y>0.1)&&(x-y<0.6))
{
	double a = (x-y-0.1);
	double b = (x-y-0.6);
	double c = 10;
	int k = 2;
	dx = (k-1)*c*(std::pow(a,k-1)*std::pow(b,k)+std::pow(a,k)*std::pow(b,k-1));
	dy = -(k-1)*c*(std::pow(a,k-1)*std::pow(b,k)+std::pow(a,k)*std::pow(b,k-1));

}*/
		else 
	{
		dx =0;
		dy =0;
	}	


};

 double CustomInitialCondition::value(double x, double y) const {
       
 		 double result = 0.0;
	/*	double radius = Hermes::sqrt((x-1)*(x-1)+y*y);
	if((radius> 0.1)&&(radius<0.9))
		{
		
		double arg = PI*(radius-0.5)/0.8;
		result= Hermes::cos(arg);
		}
		int k = 4;
		return std::pow(result,k);*/

	if((x-y>0)&&(x-y<0.5))
		{
			double arg = PI*(x-y-0.25)*2.;
			result= std::cos(arg);
			int k = 2;
			result =  std::pow(result,k);		
		}

/*
if((x-y>0.1)&&(x-y<0.6))
{
	double a = (x-y-0.1);
	double b = (x-y-0.6);
	double c = 10;
	int k = 2;
	result =std::pow(a,k)*std::pow(b,k)*c;

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

//----------bdry
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
{ double result = 0.;
  	if((x-y>0)&&(x-y<0.5))
		{
			double arg = PI*(x-y-0.25)*2.;
			result= std::cos(arg);
			int k = 2;
			result =  std::pow(result,k);		
		}
	return result;
}
