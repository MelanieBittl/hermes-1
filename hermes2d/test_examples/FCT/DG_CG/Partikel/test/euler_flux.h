#ifndef EULER_FLUX_H
#define EULER_FLUX_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


//--------------------Euler-Fluxes: Matizen A1 ( = dF_x/dU) ,A2 (= dF_y/dU)-------------------------
class EulerFluxes
{
public:
  EulerFluxes(double gamma);

  typedef double (*A_fn)(double, double, double, double, double);

  
  double A_g(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k);

double A_p(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k);

  protected:
    double gamma;
	A_fn* A_table[2];
	A_fn* A_p_table[2];
};
#endif
