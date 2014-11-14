#include "definitions.h"
const double vel_x =1.;
const double vel_y =1.;
double calc_abs_v(Element* e)
{
Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_TRIANGLE;
if(e->is_quad()) mode = HERMES_MODE_QUAD;
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
					double v_x =vel_x;
					double v_y =vel_y;
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}


CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_exact, bool DG, bool surf) : WeakForm<double>(1)
{
	
		this->set_ext(sln_exact);
   if(surf) {	
		add_matrix_form_surf(new CustomMatrixFormSurface(0, 0)); 
	add_vector_form_surf(new CustomVectorFormSurface(0) );
	}
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
  for (int i = 0; i < n; i++)
  {
	Real v_x = Real(vel_x);
	Real v_y = Real(vel_y);
   Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
//outlet
   result -= wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Scalar(0), a_dot_n) * v->val[i];

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
  for (int i = 0; i < n; i++) 
  {
	Real v_x = Real(vel_x);
	Real v_y = Real(vel_y);
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
      if(u->fn_central == NULL)
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(Real(0), u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Real(0), a_dot_n) * jump_v;
       
      
  }
  return -result;
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
  for (int i = 0; i < n; i++) {
	double v_x =vel_x;
 	double v_y =vel_y;
    double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
if(e->y[i]==0)
	 result -= wt[i] * exact->val[i] * v->val[i] * a_dot_n;
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
  { 			Real v_x = Real(vel_x);
	Real v_y = Real(vel_y);
  	result += wt[i] * (u->val[i] *(v->dx[i] * v_x + v->dy[i] * v_y ));
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
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem); //Hermes::sqrt(v_x*v_x+v_y*v_y);

   double result = double(0);
  for (int i = 0; i < n; i++)
	{		double v_x =vel_x;
 	double v_y =vel_y;
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
	double v_x =vel_x;
 	double v_y =vel_y;
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
	double v_x =vel_x;
 	double v_y =vel_y;
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * a_dot_n;
		}
  return result;
}


//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{		dx =0;
		dy =0;
		
		int k =2;
				double arg = PI*(x-y-0.4)/0.8;
				//double arg= PI*(x-0.4)/0.8;
  if(all){  
	if((x-y>0)&&(x-y<0.8))
	//if((x>0.)&&(x<0.8))
	{	dx= -std::sin(arg)*PI/0.8;
		dy=	std::sin(arg)*PI/0.8;		
		dx*= k*std::pow(std::cos(arg),k-1);
		dy *=k*std::pow(std::cos(arg),k-1);
	}
  }else{	  
		if(((x-y>0)&&(x-y<0.8))&&(x<0.8))
		//if(((x>0.)&&(x<0.8))&&(y<0.8))
		{
			dx= -std::sin(arg)*PI/0.8;
			dy=	std::sin(arg)*PI/0.8;
			dx*= k*std::pow(std::cos(arg),k-1);
			dy *=k*std::pow(std::cos(arg),k-1);
		}	  
  }
};

 double CustomInitialCondition::value(double x, double y) const 
{
      double result = 0.0; 	int k =2;
						double arg = PI*(x-y-0.4)/0.8;
				//double arg= PI*(x-0.4)/0.8;
		if(all){
	if((x-y>0)&&(x-y<0.8))
	//if((x>0.)&&(x<0.8))
			{	result= std::cos(arg);			
				result= std::pow(result,k);		
			}
		}else{			
		if(((x-y>0)&&(x-y<0.8))&&(x<0.8))
		//if(((x>0.)&&(x<0.8))&&(y<0.8))
			{
				result= std::cos(arg);
				result= std::pow(result,k);		
			}
			
		}
	
return result;

};

 Ord CustomInitialCondition::ord(double x, double y)  const 
 {
      return Ord(10);
	};
 MeshFunction<double>* CustomInitialCondition::clone() const
	{
		return new CustomInitialCondition(this->mesh,this->all);

	}



