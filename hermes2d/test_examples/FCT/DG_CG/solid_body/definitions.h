#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "hermes2d.h"
 #define PI (3.141592653589793) 

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(double time_step, double theta,  MeshFunctionSharedPtr<double> sln_prev_time, std::string inlet, MeshSharedPtr mesh, bool all = false, bool DG = true);
  WeakForm<double>* clone() const;

	
private:
  class CustomMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVol(int i, int j, double time_step, double theta) : MatrixFormVol<double>(i, j), time_step(time_step), theta(theta) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;
    
    double time_step;
    double theta;
  };

  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i, double time_step, double theta) : VectorFormVol<double>(i), time_step(time_step), theta(theta) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;

    double time_step;
    double theta;
  };

  class CustomMatrixFormSurface : public MatrixFormSurf<double>
  {
  public:
    CustomMatrixFormSurface(int i, int j) : MatrixFormSurf<double>(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormSurf<double>* clone() const;

  };

  class CustomMatrixFormInterface : public MatrixFormDG<double>
  {
  public:
    CustomMatrixFormInterface(int i, int j) : MatrixFormDG<double>(i, j) 
    {
    };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, DiscontinuousFunc<Scalar> **u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const;

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;

  };

  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface(int i, std::string inlet) : VectorFormSurf<double>(i), inlet(inlet) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;

    template<typename Real>
    Real F(Real x, Real y) const;

    template<typename Real, typename Scalar>
    Scalar g(std::string ess_bdy_marker) const;
    
    // Member.
    std::string inlet;
  };
  
  double calculate_a_dot_v(double x, double y, double vx, double vy) const;

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const;

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;

  MeshSharedPtr mesh;
};







//---------------mass-matrix/tau-----------

class CustomMatrixFormVolMassmatrix : public MatrixFormVol<double>   

{
  public:
    CustomMatrixFormVolMassmatrix(int i, int j, double time_step) 
      : MatrixFormVol<double>(i, j), time_step(time_step) { };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const ;

    // Members.  
    double time_step;
};





class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix(double time_step);
	
};

//---------------Convection-----------
class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolConvection(int i, int j) 
    : MatrixFormVol<double>(i, j) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, Func<Ord> **ext) const;  

    MatrixFormVol<double>* clone() const;

};



class CustomWeakFormConvection : public WeakForm<double>    //Konvektion
{
public:
  CustomWeakFormConvection();
	  
};


//--------------------Error_calculation--------------

class CustomNormFormVol : public NormFormVol<double>
{
public:
  CustomNormFormVol(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class StreamlineDiffusionNorm : public NormFormVol<double>
{
public:
  StreamlineDiffusionNorm(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class CustomNormFormSurf : public NormFormSurf<double>
{
public:
  CustomNormFormSurf(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class CustomNormFormDG : public NormFormDG<double>
{
public:
  CustomNormFormDG(int i, int j);

  virtual double value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const;
};



//------------------- Initial condition ----------------

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(Ord x, Ord y) const ;

   MeshFunction<double>* clone() const;
};




#endif

