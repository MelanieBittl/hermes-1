#ifndef EULER_FLUX_H
#define EULER_FLUX_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


//--------------------Euler-Fluxes: Matizen A1 ( = dF_x/dU) ,A2 (= dF_y/dU)-------------------------
class EulerFluxes
{
public:
  EulerFluxes(double kappa) : kappa(kappa) {}

  template<typename Scalar>
  Scalar A_1_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {

    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_1_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_2_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (kappa - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3. - kappa) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa - 1.);
  }

  template<typename Scalar>
  Scalar A_2_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_2_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (kappa - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3.0 - kappa) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa - 1.);
  }

  template<typename Scalar>
  Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((rho_v_x / rho) * (((kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (kappa * energy / rho)));
  }

  template<typename Scalar>
  Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((kappa * energy / rho) - (kappa - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (kappa - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    /*return Scalar(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));*/

		return Scalar(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho))));


  }

  template<typename Scalar>
  Scalar A_2_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((energy / rho) + (1 / rho) * (kappa - 1.0) * ( energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - kappa) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(kappa * rho_v_y / rho);
  }
  protected:
    double kappa;
};
#endif
