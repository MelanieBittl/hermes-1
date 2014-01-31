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
  CustomInitialCondition_rho(MeshSharedPtr mesh, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma)  {};
   ~CustomInitialCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
};

class CustomInitialCondition_rho_v_x : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho_v_x(MeshSharedPtr mesh, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma)  {};
   ~CustomInitialCondition_rho_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
};



class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(MeshSharedPtr mesh, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma) {};
   ~CustomInitialCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double gamma;
};
//--------Boundary Condition---------
class BoundaryCondition_rho : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho(MeshSharedPtr mesh,double time, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma), time(time)  {};
   ~BoundaryCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	void set_time(double t){time = t;}
	double gamma;
	double time;
};
class BoundaryCondition_rho_e : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho_e(MeshSharedPtr mesh,double time, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma), time(time)  {};
   ~BoundaryCondition_rho_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	void set_time(double t){time = t;}
	double gamma;
	double time;
};
class BoundaryCondition_rho_v_x : public ExactSolutionScalar<double>
{
public:
  BoundaryCondition_rho_v_x(MeshSharedPtr mesh,double time, double gamma) : ExactSolutionScalar<double>(mesh), gamma(gamma), time(time)  {};
   ~BoundaryCondition_rho_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	void set_time(double t){time = t;}
	double gamma;
	double time;
};


#endif
