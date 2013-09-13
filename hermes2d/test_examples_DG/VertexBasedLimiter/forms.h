#ifndef __VERTEX_BASED_FORMS
#define __VERTEX_BASED_FORMS

typedef double (*scalar_product_with_advection_direction)(double x, double y, double vx, double vy);
extern scalar_product_with_advection_direction advection_term;

#pragma region Time derivative forms

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

#pragma endregion

#pragma region Convection forms

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

class CustomVectorFormVolConvection : public VectorFormVol<double>   
{
public:

  CustomVectorFormVolConvection(int i, int ext_i, double multiplier = 1.) : VectorFormVol<double>(i), ext_i(ext_i), multiplier(multiplier)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * ext[ext_i]->val[i] * advection_term(e->x[i], e->y[i], v->dx[i], v->dy[i]);
    return multiplier * result * wf->get_current_time_step();
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
  double multiplier;
};

class CustomMatrixFormInterfaceConvection : public MatrixFormDG<double>
{
public:
  CustomMatrixFormInterfaceConvection(int i, int j, bool on_K_in, bool on_K_out, bool local) : MatrixFormDG<double>(i, j), on_K_in(on_K_in), on_K_out(on_K_out), local(local)
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
        result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? u->val_neighbor[i] : 0.) * jump_v;
      if(a_dot_n >= 0)
        result += wt[i] * a_dot_n * ((u->fn_central == NULL) ? 0. : u->val[i]) * jump_v;
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
    return new CustomMatrixFormInterfaceConvection(*this);
  }
  bool on_K_in;
  bool on_K_out;
  bool local;
};

class CustomVectorFormInterfaceConvection : public VectorFormDG<double>
{
public:
  CustomVectorFormInterfaceConvection(int i, int ext_i, bool on_K_in, bool on_K_out, double multiplier = -1.) : VectorFormDG<double>(i), ext_i(ext_i), on_K_in(on_K_in), on_K_out(on_K_out), multiplier(multiplier)
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
    return multiplier * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return ext[ext_i]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterfaceConvection(*this);
  }

  int ext_i;
  bool on_K_in;
  bool on_K_out;
  double multiplier;
};

class CustomMatrixFormSurfConvection : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurfConvection(int i, int j)
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
    return new CustomMatrixFormSurfConvection(*this);
  }

};

class CustomVectorFormSurfConvection : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurfConvection(int i, int ext_bnd, bool on_K_in, bool on_K_out, double multiplier = -1.) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), on_K_in(on_K_in), on_K_out(on_K_out), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double a_dot_n = advection_term(e->x[i], e->y[i], e->nx[i], e->ny[i]);

      if(a_dot_n >= 0 && on_K_out)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
      if(a_dot_n < 0 && on_K_in)
        result += wt[i] * a_dot_n * ext[this->ext_bnd]->val[i] * v->val[i];
    }

    return multiplier * result * wf->get_current_time_step();
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return ext[this->ext_bnd]->val[0] * v->val[0];
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurfConvection(*this);
  }

  int ext_bnd;
  bool on_K_in;
  bool on_K_out;
  double multiplier;
};

#pragma endregion

#pragma region Diffusion forms

class CustomMatrixFormVolDiffusion : public MatrixFormVol<double>   
{
public:
  CustomMatrixFormVolDiffusion(int i, int j, double diffusivity) : MatrixFormVol<double>(i, j), diffusivity(diffusivity) {}


  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0];
  }

  MatrixFormVol<double>* clone() const
  {
    return new CustomMatrixFormVolDiffusion(*this);
  }

  double diffusivity;
};

class CustomVectorFormVolDiffusion : public VectorFormVol<double>   
{
public:

  CustomVectorFormVolDiffusion(int i, int ext_i, double diffusivity, double multiplier = -1.) : VectorFormVol<double>(i), ext_i(ext_i), diffusivity(diffusivity), multiplier(multiplier)
  {
  }

  double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
  {
    double result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext[ext_i]->dx[i] * v->dx[i] + ext[ext_i]->dy[i] * v->dy[i]);
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const 
  {
    return ext[ext_i]->val[0] * v->dx[0] * e->x[0];
  }

  VectorFormVol<double>* clone() const
  {
    return new CustomVectorFormVolDiffusion(*this);
  }
  int ext_i;
  double diffusivity;
  double multiplier;
};

class CustomMatrixFormInterfaceDiffusion : public MatrixFormDG<double>
{
public:
  CustomMatrixFormInterfaceDiffusion(int i, int j, bool local, double diffusivity, double s, double sigma) : MatrixFormDG<double>(i, j), local(local), diffusivity(diffusivity), s(s), sigma(sigma)
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
    double* dx = (u->fn_central == NULL) ? u->dx_neighbor : u->dx;
    double* dy = (u->fn_central == NULL) ? u->dy_neighbor : u->dy;
    double* dx_test = (v->fn_central == NULL) ? v->dx_neighbor : v->dx;
    double* dy_test = (v->fn_central == NULL) ? v->dy_neighbor : v->dy;
    for (int i = 0; i < n; i++) 
    {
      double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
      result -= 0.5 * wt[i] * (dx[i] * e->nx[i] + dy[i] * e->ny[i]) * jump_v;
      double jump_u = (u->fn_central == NULL ? -u->val_neighbor[i] : u->val[i]);
      result -= 0.5 * wt[i] * (dx_test[i] * e->nx[i] + dy_test[i] * e->ny[i]) * jump_u * s;

      result += wt[i] * jump_v * jump_u * sigma;
    }
    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
  }

  MatrixFormDG<double>* clone() const
  {
    return new CustomMatrixFormInterfaceDiffusion(*this);
  }
  bool local;
  double diffusivity;
  double s;
  double sigma;
};

class CustomVectorFormInterfaceDiffusion : public VectorFormDG<double>
{
public:
  CustomVectorFormInterfaceDiffusion(int i, int ext_i, double diffusivity, double s, double sigma, double multiplier = -1.) : VectorFormDG<double>(i), ext_i(ext_i), diffusivity(diffusivity), s(s), sigma(sigma), multiplier(multiplier)
  {
  };

  double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
    Geom<double> *e, DiscontinuousFunc<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = .5 * (ext[this->ext_i]->dx[i] + ext[this->ext_i]->dx_neighbor[i]);
      double dy = .5 * (ext[this->ext_i]->dy[i] + ext[this->ext_i]->dy_neighbor[i]);
      result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      double dx_test = .5 * v->dx[i];
      double dy_test = .5 * v->dy[i];
      double jump = ext[this->ext_i]->val[i] - ext[this->ext_i]->val_neighbor[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * jump * s;
      result += wt[i] * sigma * jump * v->val[i];
    }
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
  {
    return ext[this->ext_i]->val[0] * v->val[0];
  }

  VectorFormDG<double>* clone() const
  {
    return new CustomVectorFormInterfaceDiffusion(*this);
  }

  int ext_i;
  double s, sigma;
  double diffusivity;
  double multiplier;
};

class CustomMatrixFormSurfDiffusion : public MatrixFormSurf<double>
{
public:
  CustomMatrixFormSurfDiffusion(int i, int j, double diffusivity, double s, double sigma, std::string inlet)
    : MatrixFormSurf<double>(i, j), diffusivity(diffusivity), s(s), sigma(sigma)
  {
    this->set_area(inlet);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double>* u, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = u->dx[i];
      double dy = u->dy[i];
      result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      double dx_test = v->dx[i];
      double dy_test = v->dy[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * u->val[i] * s;
      result += wt[i] * sigma * v->val[i] * u->val[i];
    }
    return result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord>* u, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return u->val[0] * v->val[0] * e->x[0];
  }

  MatrixFormSurf<double>* clone() const
  {
    return new CustomMatrixFormSurfDiffusion(*this);
  }
  double diffusivity;
  double s, sigma;
};

class CustomVectorFormSurfDiffusion : public VectorFormSurf<double>
{
public:
  CustomVectorFormSurfDiffusion(int i, int ext_bnd, double diffusivity, double s, double sigma, std::string inlet, bool add_grad_u = true, double multiplier = -1.) : 
    VectorFormSurf<double>(i), ext_bnd(ext_bnd), diffusivity(diffusivity), s(s), sigma(sigma), add_grad_u(add_grad_u), multiplier(multiplier)
  {
    this->set_area(inlet);
  };

  double value(int n, double *wt, Func<double> **u_ext, Func<double> *v,
    Geom<double> *e, Func<double> **ext) const
  {
    double result = 0.;
    for (int i = 0; i < n; i++) 
    {
      double dx = ext[this->ext_bnd]->dx[i];
      double dy = ext[this->ext_bnd]->dy[i];
      double dx_test = v->dx[i];
      double dy_test = v->dy[i];
      if(add_grad_u)
        result -= wt[i] * (dx * e->nx[i] + dy * e->ny[i]) * v->val[i];
      result -= wt[i] * (dx_test * e->nx[i] + dy_test * e->ny[i]) * ext[this->ext_bnd]->val[i] * s;
      result += wt[i] * sigma * v->val[i] * ext[this->ext_bnd]->val[i];
    }
    return multiplier * result * wf->get_current_time_step() * diffusivity;
  }

  Ord ord(int n, double *wt, Func<Ord> **u_ext, Func<Ord> *v,
    Geom<Ord> *e, Func<Ord> **ext) const
  {
    return ext[this->ext_bnd]->val[0] * v->val[0];
  }

  VectorFormSurf<double>* clone() const
  {
    return new CustomVectorFormSurfDiffusion(*this);
  }

  int ext_bnd;
  double diffusivity;
  double s, sigma, multiplier;
  bool add_grad_u;
};
#pragma endregion
#endif
