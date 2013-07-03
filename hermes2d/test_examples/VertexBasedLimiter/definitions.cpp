#include "definitions.h"

CustomWeakForm::CustomWeakForm() : WeakForm<double>(1) 
{
  add_matrix_form(new DefaultMatrixFormVol<double>(0, 0));
  add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));
  add_vector_form(new CustomVectorFormVol(0));
}

CustomWeakForm::CustomMatrixFormVolConvection::CustomMatrixFormVolConvection(int i, int j) : MatrixFormVol<double>(i, j)
{
}

double CustomWeakForm::CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                            Func<double> *v, Geom<double> *e, Func<double>  **ext) const                  
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (v->dx[i] + v->dy[i]);
  return -result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                       Geom<Ord> *e, Func<Ord> **ext) const 
{
  return v->val[0] * (u->dx[0] + u->dy[0]);
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}

CustomWeakForm::CustomVectorFormVol::CustomVectorFormVol(int i) : VectorFormVol<double>(i)
{
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i] * ext[0]->val[i];
  return result;
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  return v->val[0] * ext[0]->val[0];
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVol::clone() const
{
  return new CustomWeakForm::CustomVectorFormVol(*this);
}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v,
                                                        Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  double result = 0.;
  for (int i = 0; i < n; i++) 
  {
    double a_dot_n = scalar_product_with_advection_direction(e->nx[i], e->ny[i]);

    double jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * upwind_flux(0., u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * upwind_flux(u->val[i], 0., a_dot_n) * jump_v;
  }
  return result * wf->get_current_time_step();
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v,
                                                   Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return u->val[0] * u->val_neighbor[0] * v->val[0] * v->val_neighbor[0];
}

MatrixFormDG<double>* CustomWeakForm::CustomMatrixFormInterface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormInterface(*this);
}

double CustomWeakForm::CustomMatrixFormInterface::scalar_product_with_advection_direction(double vx, double vy) const
{
  return vx + vy;
}

double CustomWeakForm::CustomMatrixFormInterface::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::CustomMatrixFormInterface::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}



CustomInitialCondition::CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
{
}

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.;
  dy = 0.;
}

double CustomInitialCondition::value(double x, double y) const 
{
  if(x < 0. && y < 0.0 && x > -1. && y > -1.)
    return 1.0;
  else
    return 0.0;
}

Ord CustomInitialCondition::ord(double x, double y) const 
{
  return Ord(1);
}

MeshFunction<double>* CustomInitialCondition::clone() const
{
  return new CustomInitialCondition(this->mesh);

}

