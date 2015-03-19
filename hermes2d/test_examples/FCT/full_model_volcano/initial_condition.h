#ifndef __INITIAL_H
#define __INITIAL_H

#include "hermes2d.h"
#include "euler_util.h"


 #define PI (3.141592653589793)   
using namespace Hermes;
using namespace Hermes::Hermes2D;


const double y_max =1.;
const double c_infty = 350;
const double L = 200;
const double T_infty = 288.15;

//double mach_volcano = 0.3;
const double v_y_g = 80/c_infty;
const double v_y_p = 80/c_infty;
const double mach_left = 0.;
const double mach_init = 0.;

const double pressure_init = 101325/(c_infty*c_infty); //std::pow(10,5); //1MPa= 10^6, 0.101325 MPa = atmospheric pressure  1 bar=0.1MPa
const double pressure_coeff= 1.;
const double temp = 1200/(T_infty);
const double temp_init = 288.15/(T_infty);
const double R_gas=287.*T_infty/(c_infty*c_infty);
const double rho_p = 2300.; 
const double rho_g = pressure_init*pressure_coeff/(temp*R_gas); //p = rho*R*T
const double alpha = 0.01; //particle ash
const double alpha_2 = 1e-6; //particle else
const double vel_coeff = 1.;



// GAMMA.
const double GAMMA = 1.4; 

// Penalty Parameter.
const double SIGMA = std::pow(10,3);

const double diameter = 1e-5/L;		//10mu m
const double c_vg = 717.5*T_infty/(c_infty*c_infty);// Pelanti Diss: (gamma-1)=R/c_v
const double c_vp = 1300.*T_infty/(c_infty*c_infty);// Pelanti Diss
const double c_pg = 1004.5*T_infty/(c_infty*c_infty);// Pelanti Diss gamma = c_p/c_v
const double mu = 1e-5/(c_infty*L);// Pelanti Diss
const double kappa_g = 0.05*T_infty/(c_infty*c_infty*c_infty*L);;// Pelanti Diss
const double Pr =c_pg*mu/kappa_g;		// Pr= cp mu/kappa_g



const double g = -9.80665*L/(c_infty*c_infty);// m/s2






//------------------- Initial condition ----------------

class CustomInitialCondition_rho : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho(MeshSharedPtr mesh, double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle)  {};
   ~CustomInitialCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
	bool particle;
};

class CustomInitialCondition_rho_v_x : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho_v_x(MeshSharedPtr mesh, double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma), particle(particle)  {};
   ~CustomInitialCondition_rho_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
bool particle;
};

class CustomInitialCondition_rho_v_y : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho_v_y(MeshSharedPtr mesh, double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma), particle(particle)  {};
   ~CustomInitialCondition_rho_v_y(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
bool particle;
};

class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(MeshSharedPtr mesh, double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle) {};
   ~CustomInitialCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
bool particle;
};
//--------Boundary Condition---------
class BoundaryCondition_rho : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho(MeshSharedPtr mesh,double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle) {};
   ~BoundaryCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;

	double gamma;
bool particle;
};
class BoundaryCondition_rho_e : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho_e(MeshSharedPtr mesh,double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle)  {};
   ~BoundaryCondition_rho_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;

	double gamma;
	bool particle;
};
class BoundaryCondition_rho_v_x : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho_v_x(MeshSharedPtr mesh,double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle) {};
   ~BoundaryCondition_rho_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	
	double gamma;
bool particle;
	
};

class BoundaryCondition_rho_v_y : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho_v_y(MeshSharedPtr mesh,double gamma, bool particle=false) : ExactSolutionScalar<double>(mesh), gamma(gamma),particle(particle) {};
   ~BoundaryCondition_rho_v_y(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	
	double gamma;
bool particle;
	
};


#endif
