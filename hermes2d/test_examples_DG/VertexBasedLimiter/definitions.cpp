#include "definitions.h"

static double upwind_flux(double u_cent, double u_neib, double a_dot_n)
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

static Ord upwind_flux(Ord u_cent, Ord u_neib, double a_dot_n)
{
  return a_dot_n * (u_cent + u_neib);
}

static double advection_term_cube(double x, double y, double vx, double vy)
{
  return vx + vy;
}

static double advection_term_solid_body_rotation(double x, double y, double vx, double vy)
{
  return vx * (0.5 - y) + vy * (x - 0.5);
}

static double advection_term_circular_convection(double x, double y, double vx, double vy)
{
  return (vx * y) + (vy * (1 - x));
}

static scalar_product_with_advection_direction advection_term;

ImplicitWeakForm::ImplicitWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet, std::string outlet) : WeakForm<double>(1)
{
  switch(solvedExample)
  {
  case AdvectedCube:
    advection_term = advection_term_cube;
    break;
  case SolidBodyRotation:
    advection_term = advection_term_solid_body_rotation;
    break;
  case CircularConvection:
    advection_term = advection_term_circular_convection;
    break;
  }

  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));

  if(add_inlet)
  {
    CustomMatrixFormSurf* mform = new CustomMatrixFormSurf(0, 0);
    mform->set_area(outlet);
    this->add_matrix_form_surf(mform);
  }

  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));

  if(add_inlet)
  {
    CustomVectorFormSurf* vform = new CustomVectorFormSurf(0, 1);
    vform->set_area(inlet);
    this->add_vector_form_surf(vform);
  }

  add_vector_form(new CustomVectorFormVol(0, 0, 1.));
}

ExplicitWeakForm::ExplicitWeakForm(SolvedExample solvedExample, TimeSteppingType timeSteppingType, int explicitSchemeStep, bool add_inlet, std::string inlet, std::string outlet) : WeakForm<double>(1) 
{
  switch(solvedExample)
  {
  case AdvectedCube:
    advection_term = advection_term_cube;
    break;
  case SolidBodyRotation:
    advection_term = advection_term_solid_body_rotation;
    break;
  case CircularConvection:
    advection_term = advection_term_circular_convection;
    break;
  }

  switch(timeSteppingType)
  {
  case ExplicitEuler:
    {
      add_matrix_form(new CustomMatrixFormVol(0, 0, 1.));

      add_vector_form_DG(new CustomVectorFormInterface(0, 0));

      if(add_inlet)
      {
        CustomVectorFormSurf* bnd_form = new CustomVectorFormSurf(0, 1);
        bnd_form->set_area(inlet);
        this->add_vector_form_surf(bnd_form);

        CustomVectorFormSurfExt* ext_form = new CustomVectorFormSurfExt(0);
        ext_form->set_area(outlet);
        this->add_vector_form_surf(ext_form);
      }

      add_vector_form(new CustomVectorFormVolConvection(0));

      add_vector_form(new CustomVectorFormVol(0, 0, 1.));
    }
    break;
  case ExplicitRK:
    {
      add_vector_form_DG(new CustomVectorFormInterface(0, 0));

      add_vector_form(new CustomVectorFormVolConvection(0));

      if(add_inlet)
      {
        CustomVectorFormSurf* form = new CustomVectorFormSurf(0, 1);
        this->add_vector_form_surf(form);
      }

      if(explicitSchemeStep == 1)
      {
        add_matrix_form(new CustomMatrixFormVol(0, 0, 1.));
        add_vector_form(new CustomVectorFormVol(0, 0, 1.));
      }
      else if (explicitSchemeStep == 2)
      {
        add_matrix_form(new CustomMatrixFormVol(0, 0, 4.));
        add_vector_form(new CustomVectorFormVol(0, 0, 1.));
        add_vector_form(new CustomVectorFormVol(0, 1, 3.));
      }
      else if (explicitSchemeStep == 3)
      {
        add_matrix_form(new CustomMatrixFormVol(0, 0, 3./2.));
        add_vector_form(new CustomVectorFormVol(0, 0, 1.));
        add_vector_form(new CustomVectorFormVol(0, 1, 1./2.));
      }
    }
    break;
  }
}

double ImplicitWeakForm::CustomMatrixFormSurf::value(int n, double *wt, Func<double> **u_ext, Func<double>* u, Func<double> *v,
                                                     Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);
    result += wt[i] * u->val[i] * v->val[i] * a_dot_n;
  }
  return +result * wf->get_current_time_step();
}

Ord ImplicitWeakForm::CustomMatrixFormSurf::ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord>* u, Func<Ord> *v,
                                                Geom<Ord> *e, Func<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

MatrixFormSurf<double>* ImplicitWeakForm::CustomMatrixFormSurf::clone() const
{
  return new ImplicitWeakForm::CustomMatrixFormSurf(*this);
}

ImplicitWeakForm::CustomMatrixFormVolConvection::CustomMatrixFormVolConvection(int i, int j) : MatrixFormVol<double>(i, j)
{
}

double ImplicitWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                              Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return -result * wf->get_current_time_step();
}

Ord ImplicitWeakForm::CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                         Geom<Ord> *e, Func<Ord> **ext) const 
{
  return u->val[0] * v->dx[0];
}

MatrixFormVol<double>* ImplicitWeakForm::CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}

ImplicitWeakForm::CustomMatrixFormVol::CustomMatrixFormVol(int i, int j, double factor) : MatrixFormVol<double>(i, j), factor(factor)
{
}

double ImplicitWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * u->val[i];
  return result * this->factor;
}

Ord ImplicitWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * u->val[0];
}

MatrixFormVol<double>* ImplicitWeakForm::CustomMatrixFormVol::clone() const
{
  return new ImplicitWeakForm::CustomMatrixFormVol(*this);
}

double ImplicitWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
                                                          Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);

    if(u->fn_central == NULL)
      result += wt[i] * upwind_flux(0., u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * upwind_flux(u->val[i], 0., a_dot_n) * jump_v;
  }
  return result * wf->get_current_time_step();
}

Ord ImplicitWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
                                                     Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
}

MatrixFormDG<double>* ImplicitWeakForm::CustomMatrixFormInterface::clone() const
{
  return new ImplicitWeakForm::CustomMatrixFormInterface(*this);
}

ImplicitWeakForm::CustomVectorFormVol::CustomVectorFormVol(int i, int prev, double factor) : VectorFormVol<double>(i), prev(prev), factor(factor)
{
}

double ImplicitWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * ext[this->prev]->val[i];
  return result * this->factor;
}

Ord ImplicitWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * ext[this->prev]->val[0];
}

VectorFormVol<double>* ImplicitWeakForm::CustomVectorFormVol::clone() const
{
  return new ImplicitWeakForm::CustomVectorFormVol(*this);
}

double ImplicitWeakForm::CustomVectorFormSurf::value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
                                                     Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
  }

  return -result * wf->get_current_time_step();
}

Ord ImplicitWeakForm::CustomVectorFormSurf::ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
                                                Geom<Ord> *e, Func<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

VectorFormSurf<double>* ImplicitWeakForm::CustomVectorFormSurf::clone() const
{
  return new ImplicitWeakForm::CustomVectorFormSurf(*this);
}



ExplicitWeakForm::CustomMatrixFormVol::CustomMatrixFormVol(int i, int j, double factor) : MatrixFormVol<double>(i, j), factor(factor)
{
}

double ExplicitWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * u->val[i];
  return result * this->factor;
}

Ord ExplicitWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * u->val[0];
}

MatrixFormVol<double>* ExplicitWeakForm::CustomMatrixFormVol::clone() const
{
  return new ExplicitWeakForm::CustomMatrixFormVol(*this);
}

double ExplicitWeakForm::CustomVectorFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
                                                          Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    result += wt[i] * upwind_flux(ext[0]->val[i], ext[0]->val_neighbor[i], a_dot_n) * v->val[i];
  }
  return -result * wf->get_current_time_step();
}

Ord ExplicitWeakForm::CustomVectorFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
                                                     Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

VectorFormDG<double>* ExplicitWeakForm::CustomVectorFormInterface::clone() const
{
  return new ExplicitWeakForm::CustomVectorFormInterface(*this);
}

double ExplicitWeakForm::CustomVectorFormSurfExt::value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
                                                     Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    result += wt[i] * a_dot_n * ext[0]->val[i] * v->val[i];
  }

  return -result * wf->get_current_time_step();
}

Ord ExplicitWeakForm::CustomVectorFormSurfExt::ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
                                                Geom<Ord> *e, Func<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

VectorFormSurf<double>* ExplicitWeakForm::CustomVectorFormSurfExt::clone() const
{
  return new ExplicitWeakForm::CustomVectorFormSurfExt(*this);
}

double ExplicitWeakForm::CustomVectorFormSurf::value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
                                                     Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

    result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
  }

  return -result * wf->get_current_time_step();
}

Ord ExplicitWeakForm::CustomVectorFormSurf::ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
                                                Geom<Ord> *e, Func<Ord> **ext) const
{
  return ext[0]->val[0] * v->val[0];
}

VectorFormSurf<double>* ExplicitWeakForm::CustomVectorFormSurf::clone() const
{
  return new ExplicitWeakForm::CustomVectorFormSurf(*this);
}

ExplicitWeakForm::CustomVectorFormVolConvection::CustomVectorFormVolConvection(int i) : VectorFormVol<double>(i)
{
}

double ExplicitWeakForm::CustomVectorFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * ext[0]->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return result * wf->get_current_time_step();
}

Ord ExplicitWeakForm::CustomVectorFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
{
  return ext[0]->val[0] * v->dx[0];
}

VectorFormVol<double>* ExplicitWeakForm::CustomVectorFormVolConvection::clone() const
{
  return new CustomVectorFormVolConvection(*this);
}

ExplicitWeakForm::CustomVectorFormVol::CustomVectorFormVol(int i, int prev, double factor) : VectorFormVol<double>(i), prev(prev), factor(factor)
{
}

double ExplicitWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * ext[this->prev]->val[i];
  return result * this->factor;
}

Ord ExplicitWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * ext[this->prev]->val[0];
}

VectorFormVol<double>* ExplicitWeakForm::CustomVectorFormVol::clone() const
{
  return new ExplicitWeakForm::CustomVectorFormVol(*this);
}


InitialConditionAdvectedCube::InitialConditionAdvectedCube(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
{
}

void InitialConditionAdvectedCube::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;
  dy = 0.;
}

double InitialConditionAdvectedCube::value(double x, double y) const
{
  if(x < 0. && y < 0.0 && x > -1. && y > -1.)
    return 1.0;
  else
    return 0.0;
}

Ord InitialConditionAdvectedCube::ord(double x, double y) const 
{
  return Ord(1);
}

MeshFunction<double>* InitialConditionAdvectedCube::clone() const
{
  return new InitialConditionAdvectedCube(this->mesh);

}

void InitialConditionSolidBodyRotation::derivatives(double x, double y, double& dx, double& dy) const 
{

  double radius = 0.;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  {		
    dx = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
    dy = -std::sin(radius*M_PI)/4.0*(M_PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
  }
  else
  {			
    //cone
    x_0 = 0.5;
    y_0 = 0.25;
    radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
    if((radius< 1.0)&&(x!=x_0)) 
    { 	
      dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
      dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
    }
    else
    {
      dx=0.; dy=0.;
    }	
  }
};

double InitialConditionSolidBodyRotation::value(double x, double y) const 
{

  double result = 0.0;
  double radius;
  //hump
  double x_0 =0.25;
  double y_0= 0.5;	
  radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if( radius<= 1.0) 
  { 
    result = (1.0+ std::cos(M_PI*radius))/4.0;
    return result;	
  }
  //slotted cylinder
  x_0 = 0.5;
  y_0 = 0.75;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius <= 1) 
  { 	
    if(fabs((x-x_0))>= 0.025) return 1.0;
    if(y>=0.85) return 1.0;
  }	
  //cone
  x_0 = 0.5;
  y_0 = 0.25;
  radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
  if(radius<= 1.0) 
  { 	
    result = 1.0-radius;
  }	
  return result;
};

Ord InitialConditionSolidBodyRotation::ord(double x, double y) const 
{
  return Ord(10);
};
MeshFunction<double>* InitialConditionSolidBodyRotation::clone() const
{
  return new InitialConditionSolidBodyRotation(this->mesh);
}


void InitialConditionCircularConvection::derivatives(double x, double y, double& dx, double& dy) const 
{
  double radius = std::sqrt(std::pow(x - 1, 2.) + std::pow(y, 2.));
  double radius_dx = (2 * x - 2.) / radius;
  double radius_dy = 2 * y / radius;
  if(radius >= 0.2 && radius <= 0.4)
    dx = dy = 0.;
  else if(radius >= 0.5 && radius <= 0.8)
  {
    dx = -0.25 * std::sin(M_PI * ((radius) / 0.15)) * radius_dx;
    dx = -0.25 * std::sin(M_PI * ((radius) / 0.15)) * radius_dy;
  }
  else
    dx = dy = 0.;
};

double InitialConditionCircularConvection::value(double x, double y) const 
{
  double radius = std::sqrt(std::pow(x - 1, 2.) + std::pow(y, 2.));
  if(radius >= 0.2 && radius <= 0.4)
    return 1.;
  else if(radius >= 0.5 && radius <= 0.8)
    return 0.25 * (1. + std::cos(M_PI * ((radius - 0.65) / 0.15)));
  else
    return 0.;
};

Ord InitialConditionCircularConvection::ord(double x, double y) const 
{
  return Ord(10);
};
MeshFunction<double>* InitialConditionCircularConvection::clone() const
{
  return new InitialConditionCircularConvection(this->mesh);
}