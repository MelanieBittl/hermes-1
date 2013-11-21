#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "euler_util.h"
#include "numerical_flux.h"
#include "hermes2d.h"





  class EulerEquationsVectorFormFlux : public VectorFormDG<double>
  {
  public:
    EulerEquationsVectorFormFlux(int i, double kappa,NumericalFlux* num_flux) 
      : VectorFormDG<double>(i), num_flux(num_flux)
    {
    }

    ~EulerEquationsVectorFormFlux()
    {    };

    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const ;


    VectorFormDG<double>* clone()  const
    { 
      EulerEquationsVectorFormFlux* form = new EulerEquationsVectorFormFlux(this->i, this->num_flux->kappa, this->num_flux);
      form->wf = this->wf;
      return form;
    }

			NumericalFlux* num_flux;

  };

class EulerInterface : public WeakForm<double>
{
public:

  EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,NumericalFlux* num_flux,int num_of_equations = 4);

	~EulerInterface();

    WeakForm<double>* clone() const;


  // Members.
NumericalFlux* num_flux;



};


#endif
