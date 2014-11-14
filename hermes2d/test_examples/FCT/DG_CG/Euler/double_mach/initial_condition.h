#ifndef __INITIAL_H
#define __INITIAL_H

#include "hermes2d.h"
#include "euler_util.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


//------------------- Initial condition ----------------

class CustomInitialCondition_rho : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
};

class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(MeshSharedPtr mesh, double kappa) : ExactSolutionScalar<double>(mesh), kappa(kappa) {};
   ~CustomInitialCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double kappa;
};

class CustomInitialCondition_v_x : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_v_x(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;
  virtual MeshFunction<double>* clone() const;
 virtual Ord ord(double x, double y)  const ;
};
class CustomInitialCondition_v_y : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_v_y(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition_v_y(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;
  virtual MeshFunction<double>* clone() const;
 virtual Ord ord(double x, double y)  const ;
};


//Boundary-Condition

class CustomBoundaryCondition_rho : public ExactSolutionScalar<double>
{
public:
  CustomBoundaryCondition_rho(MeshSharedPtr mesh, double time) : ExactSolutionScalar<double>(mesh), time(time) {};
   ~CustomBoundaryCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
  
  	void set_time(double new_time){ time = new_time;}
double time;
};

class CustomBoundaryCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomBoundaryCondition_e(MeshSharedPtr mesh, double kappa, double time) : ExactSolutionScalar<double>(mesh), kappa(kappa), time(time) {};
   ~CustomBoundaryCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;
  virtual MeshFunction<double>* clone() const;
	double kappa;
	  	void set_time(double new_time){ time = new_time;}
double time;
};

class CustomBoundaryCondition_v_x : public ExactSolutionScalar<double>
{
public:
  CustomBoundaryCondition_v_x(MeshSharedPtr mesh, double time) : ExactSolutionScalar<double>(mesh), time(time) {};
   ~CustomBoundaryCondition_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;
  virtual MeshFunction<double>* clone() const;
 virtual Ord ord(double x, double y)  const ;
   	void set_time(double new_time){ time = new_time;}
double time;
};

class CustomBoundaryCondition_v_y : public ExactSolutionScalar<double>
{
public:
  CustomBoundaryCondition_v_y(MeshSharedPtr mesh, double time) : ExactSolutionScalar<double>(mesh), time(time) {};
   ~CustomBoundaryCondition_v_y(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;
  virtual MeshFunction<double>* clone() const;
 virtual Ord ord(double x, double y)  const ;
   	void set_time(double new_time){ time = new_time;}
double time;
};


#endif
