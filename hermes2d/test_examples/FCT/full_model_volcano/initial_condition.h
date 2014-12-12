#ifndef __INITIAL_H
#define __INITIAL_H

#include "hermes2d.h"
#include "euler_util.h"

 #define PI (3.141592653589793)   
using namespace Hermes;
using namespace Hermes::Hermes2D;


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
