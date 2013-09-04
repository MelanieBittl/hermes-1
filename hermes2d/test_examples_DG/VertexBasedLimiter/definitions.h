#ifndef DEF_MINE
#define DEF_MINE
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::WeakFormsH1;

enum TimeSteppingType
{
  ExplicitEuler,
  ImplicitEuler,
  ExplicitRK
};

enum SolvedExample
{
  AdvectedCube,
  SolidBodyRotation,
  CircularConvection
};

extern double upwind_flux(double u_cent, double u_neib, double a_dot_n);

extern Ord upwind_flux(Ord u_cent, Ord u_neib, double a_dot_n);

typedef double (*scalar_product_with_advection_direction)(double x, double y, double vx, double vy);
extern scalar_product_with_advection_direction advection_term;

class ImplicitWeakForm : public WeakForm<double>
{
public:
  ImplicitWeakForm(SolvedExample solvedExample, bool add_inlet = false, std::string inlet = "", std::string outlet = "");
};

class ExplicitWeakForm  : public WeakForm<double>     
{
public:
  ExplicitWeakForm(SolvedExample solvedExample, TimeSteppingType timeSteppingType, int explicitSchemeStep = 1, bool add_inlet = false, std::string inlet = "", std::string outlet = "");
};

class ExplicitWeakFormLocal  : public WeakForm<double>     
{
public:
  ExplicitWeakFormLocal(SolvedExample solvedExample, bool add_inlet, std::string inlet, std::string outlet);
};


class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  CustomMatrixFormVolConvection(int i, int j) : MatrixFormVol<double>(i, j) {}


  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * u->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
    return -result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return u->val[0] * v->dx[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVolConvection(*this);
  }
};

class CustomMatrixFormVol : public MatrixFormVol<double>   
{
public:

  CustomMatrixFormVol(int i, int j, double factor) : MatrixFormVol<double>(i, j), factor(factor)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>  **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * v->val[i] * u->val[i];
    return result * this->factor;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
  {
    return v->val[0] * u->val[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVol(*this);
  }
  double factor;
};

class CustomVectorFormVol : public VectorFormVol<double>   
{
public:

  CustomVectorFormVol(int i, int prev, double factor) : VectorFormVol<double>(i), prev(prev), factor(factor)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * v->val[i] * ext[this->prev]->val[i];
    return result * this->factor;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
  {
    return v->val[0] * ext[this->prev]->val[0];
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVol(*this);
  }
  int prev;
  double factor;
};

class CustomMatrixFormInterface : public MatrixFormDG<double>
{
public:
  CustomMatrixFormInterface(int i, int j, bool on_K_in, bool on_K_out, bool local) : MatrixFormDG<double>(i, j), on_K_in(on_K_in), on_K_out(on_K_out), local(local)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    if(local)
    {
      if(u->fn_central == NULL && v->fn_central != NULL)
        return 0.;
      if(u->fn_central != NULL && v->fn_central == NULL)
        return 0.;
    }
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

      double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);

      if(a_dot_n < 0)
      {
        if(on_K_in)
         result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? u->val_neighbor[i] : 0.) * jump_v;
        if(on_K_out)
         result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? u->val_neighbor[i] : 0.) * jump_v;
      }
      if(a_dot_n > 0)
      {
        if(on_K_in)
          result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? 0. : u->val[i]) * jump_v;
        if(on_K_out)
          result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? 0. : u->val[i]) * jump_v;
      }
    }
    return result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
  }

  MatrixFormDG<double>* clone() const
  {
    return new CustomMatrixFormInterface(*this);
  }
  bool on_K_in;
  bool on_K_out;
  bool local;
};

class CustomMatrixFormSurf : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurf(int i, int j)
    : MatrixFormSurf<double>(i, j)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double>* u, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);
      if(a_dot_n > 0)
        result += wt[i] * u->val[i] * v->val[i] * a_dot_n;
    }
    return result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord>* u, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return u->val[0] * v->val[0] * e->x[0];
  }

  MatrixFormSurf<double>* clone() const
  {
    return new CustomMatrixFormSurf(*this);
  }

};

class CustomVectorFormSurf : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurf(int i, int ext_bnd, bool on_K_in, bool on_K_out) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), on_K_in(on_K_in), on_K_out(on_K_out)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

      if(a_dot_n > 0 && on_K_out)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
      if(a_dot_n < 0 && on_K_in)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
    }

    return -result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return ext[this->ext_bnd]->val[0] * v->val[0];
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurf(*this);
  }

  int ext_bnd;
  bool on_K_in;
  bool on_K_out;
};

class CustomVectorFormInterface : public VectorFormDG<double>
{
public:
  CustomVectorFormInterface(int i, int ext_i, bool on_K_in, bool on_K_out) : VectorFormDG<double>(i), ext_i(ext_i), on_K_in(on_K_in), on_K_out(on_K_out)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);
      if(a_dot_n > 0 && on_K_out)
        result += wt[i] * a_dot_n * ext[ext_i]->val[i] * v->val[i];
      if(a_dot_n < 0 && on_K_in)
        result += wt[i] * a_dot_n * ext[ext_i]->val_neighbor[i] * v->val[i];
    }
    return -result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return ext[0]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterface(*this);
  }

  int ext_i;
  bool on_K_in;
  bool on_K_out;
};

class CustomVectorFormVolConvection : public VectorFormVol<double>   
{
public:

  CustomVectorFormVolConvection(int i, int ext_i) : VectorFormVol<double>(i), ext_i(ext_i)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * ext[ext_i]->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
    return result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return ext[ext_i]->val[0] * v->dx[0] * e->x[0];
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVolConvection(*this);
  }
  int ext_i;
};




class InitialConditionAdvectedCube : public ExactSolutionScalar<double>
{
public:
  InitialConditionAdvectedCube(MeshSharedPtr mesh);

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

class InitialConditionSolidBodyRotation : public ExactSolutionScalar<double>
{
public:
  InitialConditionSolidBodyRotation(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

class InitialConditionCircularConvection : public ExactSolutionScalar<double>
{
public:
  InitialConditionCircularConvection(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const;
};

int test();
#endif
