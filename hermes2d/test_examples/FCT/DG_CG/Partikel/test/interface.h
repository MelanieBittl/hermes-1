#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "euler_util.h"
#include "euler_flux.h"
#include "numerical_flux.h"
#include "hermes2d.h"



class EulerInterface : public WeakForm<double>
{
public:

  EulerInterface(double gamma, MeshSharedPtr mesh,MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p,NumericalFlux* num_flux,  EulerFluxes* euler_fluxes,
RiemannInvariants* riemann_invariants,int num_of_equations = 8);

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
    EulerEquationsVectorFormFlux(int entry_i, double gamma,NumericalFlux* num_flux) 
      : VectorFormDG<double>(entry_i), num_flux(num_flux), gamma(gamma), entry_i(entry_i)
    {
    }


    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const ;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const;


    VectorFormDG<double>* clone()  const
    { 
      EulerEquationsVectorFormFlux* form = new EulerInterface::EulerEquationsVectorFormFlux(this->entry_i, this->num_flux->gamma, this->num_flux);
      form->wf = this->wf;
      return form;
    }

			NumericalFlux* num_flux;
			double gamma;
			double entry_i;

  };


  class  EulerEquationsBilinearFormFlux : public MatrixFormDG<double>
  {
  public:
    EulerEquationsBilinearFormFlux(int entry_i, int entry_j, double gamma,NumericalFlux* num_flux ) : MatrixFormDG<double>(entry_i, entry_j),entry_i(entry_i), entry_j(entry_j), gamma(gamma),num_flux(num_flux)
    {
    };


    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;

double gamma;
NumericalFlux* num_flux;
	int entry_j; 
	int entry_i;


  };


};


#endif
