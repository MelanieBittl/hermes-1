#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "hermes2d.h"
 #define PI (3.141592653589793)   



using namespace Hermes;
using namespace Hermes::Hermes2D;


//---------------mass-matrix/tau-----------

class CustomMatrixFormVolMassmatrix : public MatrixFormVol<double> 
{
  public:
    CustomMatrixFormVolMassmatrix(int i, int j, double time_step) 
      : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step) { };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    // Members.  
    double time_step;
};


class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix(double time_step);
	~CustomWeakFormMassmatrix();
};

//---------------Convection-----------
class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolConvection(int i, int j) 
    : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;  

};




class CustomWeakFormConvection : public WeakForm<double>   
{
public:
  CustomWeakFormConvection();
	~CustomWeakFormConvection();  
};




//------------------- Initial condition ----------------

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(Ord x, Ord y) const ;
};



#endif

