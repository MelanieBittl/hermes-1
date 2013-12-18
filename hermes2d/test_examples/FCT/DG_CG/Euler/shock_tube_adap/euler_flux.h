#ifndef EULER_FLUX_H
#define EULER_FLUX_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

//--------------------Euler-Fluxes: Matizen A1 ( = dF_x/dU) ,A2 (= dF_y/dU)-------------------------
class EulerFluxes
{
public:
  EulerFluxes(double kappa);

  typedef double (*A_fn)(double, double, double, double, double);

  
  double A(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k);

  protected:
    double kappa;
	A_fn* A_table[2];
};
#endif
