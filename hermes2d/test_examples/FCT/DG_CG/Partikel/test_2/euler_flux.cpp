#include "euler_flux.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


//--------------------Euler-Fluxes: Matizen A1 ( = dF_x/dU) ,A2 (= dF_y/dU)-------------------------

 inline double A_1_0_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {

    return 0.;
  }

  
  inline double A_1_0_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 1.;
  }

  
  inline double A_1_0_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 0.;
  }

  
  inline double A_1_0_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 0.;
  }

  
  inline double A_2_0_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 0.;
  }

  
  inline double A_2_0_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 0.;
  }

  
  inline double A_2_0_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 1.;
  }

  
  inline double A_2_0_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return 0.;
  }

  
  inline double A_1_1_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {
    return (- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (gamma - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  inline double A_1_1_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((3. - gamma) * (rho_v_x / rho));
  }

  
  inline double A_1_1_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((1.0 - gamma) * (rho_v_y / rho));
  }

  
  inline double A_1_1_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (gamma - 1.);
  }

  
  inline double A_2_1_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (- rho_v_x * rho_v_y / (rho * rho));
  }

  
  inline double A_2_1_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (rho_v_y / rho);
  }

  
  inline double A_2_1_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (rho_v_x / rho);
  }

  
  inline double A_2_1_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return 0.;
  }

  
  inline double A_1_2_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (- rho_v_x * rho_v_y / (rho * rho));
  }

  
  inline double A_1_2_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (rho_v_y / rho);
  }

  
  inline double A_1_2_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (rho_v_x / rho);
  }

  
  inline double A_1_2_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return 0.;
  }

  
  inline double A_2_2_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (gamma - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  
  inline double A_2_2_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((1.0 - gamma) * (rho_v_x / rho));
  }

  
  inline double A_2_2_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((3.0 - gamma) * (rho_v_y / rho));
  }

  
  inline double A_2_2_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (gamma - 1.);
  }

  
  inline double A_1_3_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((rho_v_x / rho) * (((gamma - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (gamma * rho_energy / rho)));
  }

  
  inline double A_1_3_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((gamma * rho_energy / rho) - (gamma - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (gamma - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  
  inline double A_1_3_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((1.0 - gamma) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  inline double A_1_3_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (gamma * (rho_v_x / rho));
  }

  
  inline double A_2_3_0(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    /*return (- (rho_v_y * rho_energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (gamma - 1.0) 
            * (rho_energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (gamma - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));*/

		return (- (rho_v_y * rho_energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (gamma - 1.0) 
            * (rho_energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho))));


  }

  
  inline double A_2_3_1(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((1.0 - gamma) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  
  inline double A_2_3_2(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return ((rho_energy / rho) + (1 / rho) * (gamma - 1.0) * ( rho_energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - gamma) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  
  inline double A_2_3_3(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (gamma * rho_v_y / rho);
  }

//------------additional for partikels

 inline double zero_value(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {

    return 0.;
  }
 inline double one_value(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma) {

    return 1.;
  }
inline double A_1_1_0_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return -Hermes::sqr(rho_v_x / rho);
  }
inline double A_1_1_1_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (2.*(rho_v_x / rho));
  }


  inline double A_1_3_0_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return -(rho_energy / rho)*(rho_v_x / rho);
  }

  
  inline double A_1_3_1_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return  (rho_energy / rho);
  }

  
  inline double A_1_3_3_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return  (rho_v_x / rho);
  }

inline double A_2_2_0_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return -Hermes::sqr(rho_v_y / rho);
  }
inline double A_2_2_2_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (2.*(rho_v_y / rho));
  }


 inline double A_2_3_0_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
		return (- (rho_v_y * rho_energy) / (rho * rho));


  }


  
  inline double A_2_3_2_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return  (rho_energy / rho);
  }

  
  inline double A_2_3_3_p(double rho, double rho_v_x, double rho_v_y, double rho_energy , double gamma){
    return (rho_v_y / rho);
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


    static EulerFluxes::A_fn A_1_p_fn[] =
    {
		zero_value, one_value, zero_value, zero_value, 
		A_1_1_0_p, A_1_1_1_p, zero_value, zero_value, 
		A_1_2_0, A_1_2_1, A_1_2_2, zero_value,
		A_1_3_0_p, A_1_3_1_p, zero_value, A_1_3_3_p
    };

    static EulerFluxes::A_fn A_2_p_fn[] =
    {
		zero_value,zero_value, one_value, zero_value, 
		A_2_1_0, A_2_1_1, A_2_1_2, zero_value, 
		A_2_2_0_p, zero_value, A_2_2_2_p,zero_value,
		A_2_3_0_p, zero_value, A_2_3_2_p, A_2_3_3_p
    };

    static EulerFluxes::A_fn* A_p_fn_table[2] =
    {
      A_1_p_fn,
	  A_2_p_fn
    };



  
  double EulerFluxes::A_g(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k)
	{
		
		return A_table[i][j*4+k](rho, rho_v_x, rho_v_y, rho_energy, gamma);

	}


  double EulerFluxes::A_p(double rho, double rho_v_x, double rho_v_y, double rho_energy, int i, int j, int k)
	{
		
		return A_p_table[i][j*4+k](rho, rho_v_x, rho_v_y, rho_energy, gamma);

	}

EulerFluxes::EulerFluxes(double gamma) : gamma(gamma) {
	A_table[0] = A_fn_table[0];
	A_table[1] = A_fn_table[1];
	A_p_table[0] = A_p_fn_table[0];
	A_p_table[1] = A_p_fn_table[1];

}

