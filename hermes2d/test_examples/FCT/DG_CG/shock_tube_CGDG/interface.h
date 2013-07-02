#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "euler_util.h"
#include "numerical_flux.h"
#include "hermes2d.h"



 class EulerInterfaceBilinearForm: public MatrixFormDG<double>
  {
  public:
    EulerInterfaceBilinearForm(int i, int j, double kappa, EulerFluxes* fluxes, bool* cacheReady, double** P_plus_cache, double** P_minus_cache) 
      : MatrixFormDG<double>(i, j), num_flux(new StegerWarmingNumericalFlux(kappa)), cacheReady(cacheReady), P_plus_cache(P_plus_cache), P_minus_cache(P_minus_cache), fluxes(fluxes) 
    {
    }

    ~EulerInterfaceBilinearForm() 
    {
      delete num_flux;
    }
    double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, 
      DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const ;

    MatrixFormDG<double>* clone()  const;


    bool* cacheReady;
    double** P_plus_cache;
    double** P_minus_cache;
    StegerWarmingNumericalFlux* num_flux;
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
	RiemannInvariants* riemann_invariants;

  bool cacheReadyDG;
  double** P_plus_cache_DG;
  double** P_minus_cache_DG;

};


#endif
