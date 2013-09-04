#include "definitions.h"

double upwind_flux(double u_cent, double u_neib, double a_dot_n)
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord upwind_flux(Ord u_cent, Ord u_neib, double a_dot_n)
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

scalar_product_with_advection_direction advection_term;

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

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterface(0, 0, false, true, false));
  // Numerical flux - inner-edge inlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterface(0, 0, false, false));
  // Numerical flux - inner-edge inlet+outlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterface(0, 1, true, true));

  // No convection - test functions have zero gradient
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));

  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurf(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurf(0, 2, true, false));
    // Numerical flux - boundary outlet - derivatives - rhs
    this->add_vector_form_surf(new CustomVectorFormSurf(0, 1, false, true));
  }

  // Mass matrix - rhs
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

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
  
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterface(0, 0, false, true, true));
  // Numerical flux - inner-edge inlet+outlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterface(0, 0, true, true));
  // Numerical flux - inner-edge inlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterface(0, 1, true, false));

  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurf(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurf(0, 2, true, false));
    // Numerical flux - boundary outlet - mean values - rhs
    this->add_vector_form_surf(new CustomVectorFormSurf(0, 0, false, true));
  }

  // Convective term - matrix
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  // Convective term - rhs
  add_vector_form(new CustomVectorFormVolConvection(0, 0));

  // Mass matrix - rhs
  add_vector_form(new CustomVectorFormVol(0, 1, 1.));
}

/*
ExplicitWeakFormLocal::ExplicitWeakFormLocal(SolvedExample solvedExample, bool add_inlet, std::string inlet, std::string outlet)
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
  this->add_matrix_form(new ExplicitWeakFormLocal::CustomMatrixFormVolConvection(0, 0));
  this->add_matrix_form_DG(new ExplicitWeakFormLocal::CustomMatrixFormOutgoing(0, 0));
  this->add_vector_form_DG(new ExplicitWeakFormLocal::CustomVectorFormIncoming(0));

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
}

*/
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
  return Ord(20);
};
MeshFunction<double>* InitialConditionCircularConvection::clone() const
{
  return new InitialConditionCircularConvection(this->mesh);
}