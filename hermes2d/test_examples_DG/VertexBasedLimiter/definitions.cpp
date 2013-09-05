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

static double advection_term_moving_peak(double x, double y, double vx, double vy)
{
  return -(vx * y) + (vy * x);
}

static double advection_term_benchmark(double x, double y, double vx, double vy)
{
  return vx + vy;
}
scalar_product_with_advection_direction advection_term;

ImplicitWeakForm::ImplicitWeakForm(SolvedExample solvedExample, bool add_inlet, std::string inlet, std::string outlet, double diffusivity) : WeakForm<double>(1)
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
  case MovingPeak:
    advection_term = advection_term_moving_peak;
    break;
  case Benchmark:
    advection_term = advection_term_benchmark;
    break;
  }

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
  // Mass matrix - rhs
  add_vector_form(new CustomVectorFormVol(0, 0, 1.));
  add_vector_form(new CustomVectorFormVol(0, 1, 1.));

  // Convection.
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, false, true, false));
  // Numerical flux - inner-edge inlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, false, false));
  // Numerical flux - inner-edge inlet+outlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 1, true, true));
  // No convection - when test functions have zero gradient
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 2, true, false));
    // Numerical flux - boundary outlet - derivatives - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 1, false, true));
  }

  // Diffusion.
  // No Diffusion - when test functions have zero gradient
  add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
  add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, false, diffusivity, -1, 1. * 100.));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, -1, 1. * 100., inlet));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 2, diffusivity, -1, 1. * 100., inlet));
  }
}

ExplicitWeakForm::ExplicitWeakForm(SolvedExample solvedExample, TimeSteppingType timeSteppingType, int explicitSchemeStep, bool add_inlet, std::string inlet, std::string outlet, double diffusivity) : WeakForm<double>(1) 
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
  case MovingPeak:
    advection_term = advection_term_moving_peak;
    break;
  case Benchmark:
    advection_term = advection_term_benchmark;
    break;
  }

  // Mass matrix
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));

  // Convection.
  // Convective term - matrix
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  // Convective term - rhs
  add_vector_form(new CustomVectorFormVolConvection(0, 0));
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceConvection(0, 0, false, true, true));
  // Numerical flux - inner-edge inlet+outlet - mean values - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 0, true, true));
  // Numerical flux - inner-edge inlet - derivatives - rhs
  add_vector_form_DG(new CustomVectorFormInterfaceConvection(0, 1, true, false));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfConvection(0, 0));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 2, true, false));
    // Numerical flux - boundary outlet - mean values - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfConvection(0, 0, false, true));
  }

  // Diffusion.
  // Diffusive term - matrix
  add_matrix_form(new CustomMatrixFormVolDiffusion(0, 0, diffusivity));
  // Diffusive term - rhs - for mean values - zero
  add_vector_form(new CustomVectorFormVolDiffusion(0, 0, diffusivity));
  // Numerical flux - inner-edge element outlet - matrix
  add_matrix_form_DG(new CustomMatrixFormInterfaceDiffusion(0, 0, true, diffusivity, -1, 1. * 100.));
  if(add_inlet)
  {
    // Numerical flux - boundary outlet - matrix
    this->add_matrix_form_surf(new CustomMatrixFormSurfDiffusion(0, 0, diffusivity, -1, 1. * 100., inlet));
    // Numerical flux - boundary inlet - exact solution - rhs
    this->add_vector_form_surf(new CustomVectorFormSurfDiffusion(0, 2, diffusivity, -1, 1. * 100., inlet));
  }

  // Mass matrix - rhs
  add_vector_form(new CustomVectorFormVol(0, 0, 1.));
  add_vector_form(new CustomVectorFormVol(0, 1, 1.));
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


void ExactSolutionCircularConvection::derivatives(double x, double y, double& dx, double& dy) const 
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

double ExactSolutionCircularConvection::value(double x, double y) const 
{
  double radius = std::sqrt(std::pow(x - 1, 2.) + std::pow(y, 2.));
  if(radius >= 0.2 && radius <= 0.4)
    return 1.;
  else if(radius >= 0.5 && radius <= 0.8)
    return 0.25 * (1. + std::cos(M_PI * ((radius - 0.65) / 0.15)));
  else
    return 0.;
};

Ord ExactSolutionCircularConvection::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionCircularConvection::clone() const
{
  return new ExactSolutionCircularConvection(this->mesh);
}

void ExactSolutionMovingPeak::set_current_time(double time)
{
  this->current_time = time;
  double x_0 = 0.;
  double y_0 = 0.5;
  this->x_hat = x_0 * std::cos(this->current_time) - y_0 * std::sin(this->current_time);
  this->y_hat = -x_0 * std::sin(this->current_time) + y_0 * std::cos(this->current_time);
}

void ExactSolutionMovingPeak::derivatives(double x, double y, double& dx, double& dy) const 
{
  double fraction_times_pi = 1. / (4. * this->diffusivity * this->current_time);
  double fn_value = this->value(x, y);
  dx = fn_value * 2. * (x - this->x_hat) * fraction_times_pi;
  dy = fn_value * 2. * (x - this->x_hat) * fraction_times_pi;
};

double ExactSolutionMovingPeak::value(double x, double y) const 
{
  double fraction_times_pi = 1. / (4. * this->diffusivity * this->current_time);
  double r_squared = ( (x - this->x_hat) * (x - this->x_hat) ) + ( (y - this->y_hat) * (y - this->y_hat) );
  return fraction_times_pi * std::pow(2.718281828, -r_squared * fraction_times_pi) / M_PI;
};

Ord ExactSolutionMovingPeak::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionMovingPeak::clone() const
{
  return new ExactSolutionMovingPeak(this->mesh, this->diffusivity, this->current_time);
}

double ExactSolutionMovingPeak::get_current_time() const
{
  return this->current_time;
}



void InitialConditionBenchmark::derivatives(double x, double y, double& dx, double& dy) const 
{
  double x_0 = 2./15.;
  double l = 7. * std::sqrt(2.0)  / 300.;
  double fn_value = this->value(x, y);
  dx = fn_value * (-2. * (x - x_0) / l);
};

double InitialConditionBenchmark::value(double x, double y) const 
{
  double x_0 = 2./15.;
  double l = 7. * std::sqrt(2.0)  / 300.;
  return (5./7.) * std::exp(-std::pow((x - x_0) / l, 2.));
};

Ord InitialConditionBenchmark::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* InitialConditionBenchmark::clone() const
{
  return new InitialConditionBenchmark(this->mesh);
}

void ExactSolutionBenchmark::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = dy = 0.;
};

double ExactSolutionBenchmark::value(double x, double y) const 
{
  double x_0 = 2./15.;
  double l = 7. * std::sqrt(2.0)  / 300.;
  double fraction = 5./(7. * sigma());
  return fraction * std::exp(-std::pow((x - x_0 - this->time) / (l * sigma()), 2.));
};

double ExactSolutionBenchmark::sigma() const
{
  double l = 7. * std::sqrt(2.0)  / 300.;
  return std::sqrt(1. + ((4. * this->diffusivity * this->time) / (l*l)));
}

Ord ExactSolutionBenchmark::ord(double x, double y) const 
{
  return Ord(20);
};
MeshFunction<double>* ExactSolutionBenchmark::clone() const
{
  return new ExactSolutionBenchmark(this->mesh, this->diffusivity, this->time);
}

