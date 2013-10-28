#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "euler_util.h"
#include "numerical_flux.h"
#include "hermes2d.h"





  class EulerEquationsVectorFormFlux : public VectorFormDG<double>
  {
  public:
    EulerEquationsVectorFormFlux(int i, double kappa, EulerFluxes* fluxes) 
      : VectorFormDG<double>(i), num_flux(new HLLNumericalFlux(kappa)), fluxes(fluxes) 
    {
    }

    ~EulerEquationsVectorFormFlux()
    {
      delete num_flux;
    }

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double result = 0.;
      double w_L[4], w_R[4];
      for (int point_i = 0; point_i < n; point_i++)
      {
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];

        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, e->nx[point_i], e->ny[point_i]) * v->val[point_i];
      }

      return -result * wf->get_current_time_step();
    }

    VectorFormDG<double>* clone()  const
    { 
      EulerEquationsVectorFormFlux* form = new EulerEquationsVectorFormFlux(this->i, this->num_flux->kappa, this->fluxes);
      form->wf = this->wf;
      return form;
    }

    HLLNumericalFlux* num_flux;
    EulerFluxes* fluxes;
  };

class EulerInterface : public WeakForm<double>
{
public:

  EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,  int num_of_equations = 4);

	~EulerInterface();

    WeakForm<double>* clone() const;


  // Members.
  EulerFluxes* euler_fluxes;



};


#endif
