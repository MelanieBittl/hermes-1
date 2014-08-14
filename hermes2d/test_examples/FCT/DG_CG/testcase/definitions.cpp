#include "definitions.h"
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
				double v_x = 2.-y_coord[j]*y_coord[j];
				double v_y = 2- x_coord[j];	 
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}	
	delete rm;
	return Hermes::sqrt(abs_v);
}

CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,bool all, bool DG, bool SD) : WeakForm<double>(1)
{
 this->set_ext(sln_prev_time);

  if(all)
  {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
		add_vector_form(new RHS(0));
  } 
   if(DG) add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));
	 if(SD)
	{	add_matrix_form(new Streamline(0,0,mesh));
		add_vector_form(new RHS_sd(0,mesh));
		}
   
	   add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));   
		add_vector_form_surf(new CustomVectorFormSurface(0) );
}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}

//---------------convection-----------

    template<typename Real, typename Scalar>
    Scalar CustomWeakForm::CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{
			double x = e->x[i]; 
			double y = e->y[i];
				double v_x = 2.-y*y;
				double v_y = 2- x;
				double c = 1+(1+x)*(1+y)*(1+y);

		result += wt[i] * (u->val[i] *(v->val[i]*c-(v->dx[i]*v_x+ v->dy[i]*v_y)));

}
  return result;

    };

double CustomWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
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
template<typename Real, typename Scalar>
Scalar CustomWeakForm::Streamline::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);

		Real diam = e->diam;	
		//Real area = Hermes::sqrt(e->area);
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem);
for (int i = 0; i < n; i++)
{
			Real x = e->x[i]; 
			Real y = e->y[i];
				Real v_x = 2.-y*y;
				Real v_y = 2- x;	
			Real c = 1+(1+x)*(1+y)*(1+y);		
	Real tau = diam/(2.*abs_v);
  result += wt[i] *(tau*(u->dx[i] * v_x + u->dy[i] * v_y + u->val[i]*c ) 
  *(v->dx[i] * v_x + v->dy[i] * v_y ));

}


  return result;
}

double CustomWeakForm::Streamline::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::Streamline::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakForm::Streamline::clone() const
{
  return new CustomWeakForm::Streamline(*this);
}




//--------Boundary-Matrix-Form
double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
 
	for (int i = 0; i < n; i++)
	{

	//outlet
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
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
template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
			Real x = e->x[i]; 
			Real y = e->y[i];
				Real v_x = 2.-y*y;
				Real v_y = 2- x;	
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
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
//----------Dirichlet Vector surface Form
double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
 Func<double>* exact = ext[0];	
   for (int i = 0; i < n; i++)
	{ 
		 //Dirichlet-Bdry!
				{
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
					double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
					result -= wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0., exact->val[i] , a_dot_n) * v->val[i];					
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
 	 Func<double>* exact = ext[0];
  for (int i = 0; i < n; i++)
		{ 	
			double x = e->x[i]; 
			double y = e->y[i];
			double v_x = 2.-y*y;
			double v_y = 2- x;	
			double c = 1+(1+x)*(1+y)*(1+y);
			result += wt[i] * (exact->val[i]*c + v_x*exact->dx[i] + v_y*exact->dy[i]) * v->val[i];
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


//--------rhs Vector Form streamline -----------------
double  CustomWeakForm::RHS_sd::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{  
	double result = 0;
	Func<double>* exact = ext[0];
		double diam = e->diam;
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem);
   for (int i = 0; i < n; i++)
		{ 	
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
			double tau = diam/(2.*abs_v);
			double c = 1+(1+x)*(1+y)*(1+y);
			result += wt[i] * (exact->val[i]*c + v_x*exact->dx[i] + v_y*exact->dy[i]) * tau* (v->dx[i] * v_x + v->dy[i] * v_y );
		}
 return result;
}

    Ord  CustomWeakForm::RHS_sd::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
	return Ord(10);
}

    VectorFormVol<double>*  CustomWeakForm::RHS_sd::clone() const
{

	return new CustomWeakForm::RHS_sd(*this);
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
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
    result += wt[i] * (v_x*u->dx[i] + v_y*u->dy[i]) * (v_x*v->dx[i] + v_y*v->dy[i]);
	}
  return (result*diam/(abs_v));
	

}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
   double result = double(0);
  for (int i = 0; i < n; i++)
	{
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
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
			double x = e->x[i]; 
			double y = e->y[i];
							double v_x = 2.-y*y;
				double v_y = 2- x;	
		double a_dot_n = std::abs(v_x*e->nx[i]+ v_y* e->ny[i]);
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]) * a_dot_n;
		}
  return result;
}
//------------------- Initial condition ----------------

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
       double result = 0.0;

		double arg =(1+x)*(1+y)*(1+y)/8.;
		 result = std::cos(PI*arg); 
		 dx = result*(1+y)*(1+y)*PI/8.;
		 dy = result*(1+x)*(1+y)*PI/4.;


}

 double CustomInitialCondition::value(double x, double y) const 
{
       
  double result = 0.0;

double arg =(1+x)*(1+y)*(1+y)/8.;
		 result = (1.0+ std::sin(PI*arg));
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
