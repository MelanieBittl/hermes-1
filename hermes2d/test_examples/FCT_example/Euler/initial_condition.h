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
  CustomInitialCondition_rho(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(Ord x, Ord y) const ;
};

class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(Mesh* mesh, double kappa) : ExactSolutionScalar<double>(mesh), kappa(kappa) {};
   ~CustomInitialCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(Ord x, Ord y) const ;
	double kappa;
};



/// Custom boundary condition
class CustomBC_v_x : public EssentialBoundaryCondition<double> {
public:
  CustomBC_v_x(std::string marker, double kappa) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), kappa(kappa)
  {
    markers.push_back(marker);
  }

  ~CustomBC_v_x() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

	double kappa;
};

class CustomBC_v_y : public EssentialBoundaryCondition<double> {
public:
  CustomBC_v_y(std::string marker, double kappa) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), kappa(kappa)
  {
    markers.push_back(marker);
  }

  ~CustomBC_v_y() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

	double kappa;

};
class CustomBC_rho : public EssentialBoundaryCondition<double> {
public:
  CustomBC_rho(std::string marker, double kappa) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), kappa(kappa)
  {
    markers.push_back(marker);
  }

  ~CustomBC_rho() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

	double kappa;

};

class CustomBC_e : public EssentialBoundaryCondition<double> {
public:
  CustomBC_e(std::string marker, double kappa) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), kappa(kappa)
  {
    markers.push_back(marker);
  }

  ~CustomBC_e() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

	double kappa;

};


#endif
