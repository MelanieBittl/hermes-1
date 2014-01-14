#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "euler_util.h"
#include "euler_flux.h"
#include "numerical_flux.h"
#include "hermes2d.h"



class EulerInterface : public WeakForm<double>
{
public:

  EulerInterface(double kappa, MeshSharedPtr mesh,MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,NumericalFlux* num_flux,  EulerFluxes* euler_fluxes,
RiemannInvariants* riemann_invariants,int num_of_equations = 4);

	~EulerInterface();

    WeakForm<double>* clone() const;


  // Members.
NumericalFlux* num_flux;
MeshSharedPtr mesh;
  EulerFluxes* euler_fluxes;
RiemannInvariants* riemann_invariants;


protected:

  class EulerEquationsVectorFormFlux : public VectorFormDG<double>
  {
  public:
    EulerEquationsVectorFormFlux(int i, double kappa,NumericalFlux* num_flux) 
      : VectorFormDG<double>(i), num_flux(num_flux), kappa(kappa)
    {
    }


    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const ;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const;


    VectorFormDG<double>* clone()  const
    { 
      EulerEquationsVectorFormFlux* form = new EulerInterface::EulerEquationsVectorFormFlux(this->i, this->num_flux->kappa, this->num_flux);
      form->wf = this->wf;
      return form;
    }

			NumericalFlux* num_flux;
			double kappa;

  };


  class  EulerEquationsBilinearFormFlux : public MatrixFormDG<double>
  {
  public:
    EulerEquationsBilinearFormFlux(int i, int j, double kappa ) : MatrixFormDG<double>(i, j), kappa(kappa)
    {
    };


    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;

double kappa;


  };


};


#endif
