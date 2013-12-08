#include "euler_flux.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


//--------------------Euler-Fluxes: Matizen A1 ( = dF_x/dU) ,A2 (= dF_y/dU)-------------------------

 inline double A_1_0_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {

    return 0.;
  }

  
  inline double A_1_0_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 1.;
  }

  
  inline double A_1_0_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 0.;
  }

  
  inline double A_1_0_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 0.;
  }

  
  inline double A_2_0_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 0.;
  }

  
  inline double A_2_0_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 0.;
  }

  
  inline double A_2_0_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 1.;
  }

  
  inline double A_2_0_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return 0.;
  }

  
  inline double A_1_1_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa) {
    return (- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (kappa - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  inline double A_1_1_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((3. - kappa) * (rho_v_x / rho));
  }

  
  inline double A_1_1_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((1.0 - kappa) * (rho_v_y / rho));
  }

  
  inline double A_1_1_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (kappa - 1.);
  }

  
  inline double A_2_1_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (- rho_v_x * rho_v_y / (rho * rho));
  }

  
  inline double A_2_1_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (rho_v_y / rho);
  }

  
  inline double A_2_1_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (rho_v_x / rho);
  }

  
  inline double A_2_1_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return 0.;
  }

  
  inline double A_1_2_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (- rho_v_x * rho_v_y / (rho * rho));
  }

  
  inline double A_1_2_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (rho_v_y / rho);
  }

  
  inline double A_1_2_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (rho_v_x / rho);
  }

  
  inline double A_1_2_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return 0.;
  }

  
  inline double A_2_2_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (kappa - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  inline double A_2_2_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((1.0 - kappa) * (rho_v_x / rho));
  }

  
  inline double A_2_2_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((3.0 - kappa) * (rho_v_y / rho));
  }

  
  inline double A_2_2_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (kappa - 1.);
  }

  
  inline double A_1_3_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((rho_v_x / rho) * (((kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (kappa * rho_energy / rho)));
  }

  
  inline double A_1_3_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((kappa * rho_energy / rho) - (kappa - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (kappa - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  
  inline double A_1_3_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  inline double A_1_3_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (kappa * (rho_v_x / rho));
  }

  
  inline double A_2_3_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    /*return (- (rho_v_y * rho_energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (rho_energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (kappa - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));*/

		return (- (rho_v_y * rho_energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (kappa - 1.0) 
            * (rho_energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho))));


  }

  
  inline double A_2_3_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((1.0 - kappa) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  inline double A_2_3_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return ((rho_energy / rho) + (1 / rho) * (kappa - 1.0) * ( rho_energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - kappa) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  
  inline double A_2_3_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double kappa){
    return (kappa * rho_v_y / rho);
  }

///-----------------------
    static EulerFluxes::A_fn A_1_fn[] =
    {
		A_1_0_0, A_1_0_1, A_1_0_2, A_1_0_3, 
		A_1_1_0, A_1_1_1, A_1_1_2, A_1_1_3, 
		A_1_2_0, A_1_2_1, A_1_2_2, A_1_2_3,
		A_1_3_0, A_1_3_1, A_1_3_2, A_1_3_3
    };

    static EulerFluxes::A_fn A_2_fn[] =
    {
		A_2_0_0, A_2_0_1, A_2_0_2, A_2_0_3, 
		A_2_1_0, A_2_1_1, A_2_1_2, A_2_1_3, 
		A_2_2_0, A_2_2_1, A_2_2_2, A_2_2_3,
		A_2_3_0, A_2_3_1, A_2_3_2, A_2_3_3
    };

    static EulerFluxes::A_fn* A_fn_table[2] =
    {
      A_1_fn,
	  A_2_fn
    };


  
  double EulerFluxes::A(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k)
	{
		
		return A_table[i][j*4+k](rho, rho_v_x, rho_v_y, rho_energy, kappa);

	}

EulerFluxes::EulerFluxes(double kappa) : kappa(kappa) {
	A_table[0] = A_fn_table[0];
	A_table[1] = A_fn_table[1];

}

