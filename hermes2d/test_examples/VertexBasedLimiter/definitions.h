#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::WeakFormsH1;

class CustomWeakForm  : public WeakForm<double>     
{
public:
  CustomWeakForm();

  class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
  {
  public:
    CustomMatrixFormVolConvection(int i, int j);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    MatrixFormVol<double>* clone() const;
  };

  class CustomVectorFormVol : public VectorFormVol<double>   
  {
  public:
    CustomVectorFormVol(int i);

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const;  

    VectorFormVol<double>* clone() const;
  };

  class CustomMatrixFormInterface : public MatrixFormDG<double>
  {
  public:
    CustomMatrixFormInterface(int i, int j) : MatrixFormDG<double>(i, j) 
    {
    };

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;

    double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

    Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;

    double scalar_product_with_advection_direction(double vx, double vy) const;
  };
};

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh);

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};