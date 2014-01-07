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
  CustomInitialCondition_rho(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)   const ;
  virtual MeshFunction<double>* clone() const;
};

class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(MeshSharedPtr mesh, double kappa) : ExactSolutionScalar<double>(mesh), kappa(kappa) {};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)   const ;
  virtual MeshFunction<double>* clone() const;
	double kappa;
};


class CustomInitialCondition_p : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_p(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)   const ;
	double kappa;
  virtual MeshFunction<double>* clone() const;
};

#endif
