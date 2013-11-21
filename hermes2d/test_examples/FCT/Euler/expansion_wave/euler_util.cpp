#include "euler_util.h"
#include "limits.h"
#include <limits>

// Calculates energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa)
{
  double to_return = pressure/(kappa - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
 // if(std::abs(to_return) < 1E-12 || to_return < 0.0)
   // return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = (kappa - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = std::sqrt(kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy, kappa) / rho);
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

  // Calculates enthalpy.
double QuantityCalculator::enthalpy(double rho, double rho_v_x, double rho_v_y, double energy, double kappa){

double to_return= (energy+ QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y,energy,kappa))/rho;
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
 if(to_return < 0.0)
   return 0.0;

return to_return;
}


void MachNumberFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(kappa * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / values.at(0)[i]);
}


void PressureFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (kappa - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}

void VelocityFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = values.at(coord)[i]/values.at(0)[i];
}


void EntropyFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / p_ext)
      / Hermes::pow((values.at(0)[i] / rho_ext), kappa));
}


//-------------------RiemannInvariants-----------------------
double RiemannInvariants::get_w1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - 2*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa)/(kappa-1));
}

double RiemannInvariants::get_w2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return QuantityCalculator::calc_pressure(rho, rho_v_x,rho_v_y, rho_energy, kappa)/Hermes::pow(rho,kappa);
}

double RiemannInvariants::get_w3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double t_x, double t_y){
		return  (rho_v_x*t_x/rho + rho_v_y*t_y/rho);
}

double RiemannInvariants::get_w4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +2*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa)/(kappa-1));
}

bool RiemannInvariants::solid_wall(double rho, double rho_v_x, double rho_v_y, double n_x, double n_y){
	if(fabs(rho_v_x*n_x/rho + rho_v_y*n_y/rho)<1E-12) 
				return true;
	else false;

}


double RiemannInvariants::get_ev1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa));
}

double RiemannInvariants::get_ev2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa));
}

double RiemannInvariants::get_pressure(double w_1, double w_2, double w_3, double w_4){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return (c*c*rho/kappa);
}

double RiemannInvariants::get_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x){
	return (w_1+w_4)*n_x/2.+w_3*t_x;
}
double RiemannInvariants::get_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y){
	return (w_1+w_4)*n_y/2.+w_3*t_y;

}
double RiemannInvariants::get_speed_sound(double w_1, double w_4){
									return (kappa-1.)*(w_4-w_1)/4.;
}

double RiemannInvariants::get_rho(double w_1, double w_2, double w_3, double w_4){
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return  Hermes::pow(c*c/(kappa*w_2),(1/(kappa-1)));
}
double RiemannInvariants::get_rho_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double v_x = RiemannInvariants::get_v_x(w_1,w_2, w_3, w_4, n_x, t_x);
	return v_x*rho;

}
double RiemannInvariants::get_rho_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double v_y = RiemannInvariants::get_v_y(w_1,w_2, w_3, w_4, n_y, t_y);
		return v_y*rho;

}
double RiemannInvariants::get_energy(double w_1, double w_2, double w_3, double w_4, double n_x, double n_y, double t_x, double t_y){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double p = RiemannInvariants::get_pressure(w_1,w_2, w_3, w_4);
	double v_x = RiemannInvariants::get_v_x(w_1,w_2, w_3, w_4, n_x, t_x);
	double v_y = RiemannInvariants::get_v_y(w_1,w_2, w_3, w_4, n_y, t_y);
	double v_norm = v_x*v_x+v_y*v_y;	
	return (p/(kappa-1.)+rho*v_norm/2.);

}


void RiemannInvariants::get_ghost_state(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid)
{

		double w_1,w_2,w_3,w_4;	

	if( solid==true){
		ghost_state[0] = rho;
		ghost_state[3] = rho_energy;
		ghost_state[1] = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
		ghost_state[2] = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
	}else if(bdry==1){ //supersonic outlet{
		ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=		rho_v_y;
		ghost_state[3]= rho_energy;
	}else if(bdry==2){//supersonic inlet
		ghost_state[0] =  rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=		rho_v_y_ext;
		ghost_state[3]= rho_energy_ext;
	
	}
}

int RiemannInvariants::get_bdry_info(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid){

	double lambda_1 = RiemannInvariants::get_ev1(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_2_3 = RiemannInvariants::get_ev2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_4 = RiemannInvariants::get_ev4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		double w_1,w_2,w_3,w_4;	

	if( solid==true){
		ghost_state[0] = rho;
		ghost_state[3] = rho_energy;
		ghost_state[1] = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
		ghost_state[2] = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
			return 0;
	}else if((lambda_1>=0)&&(lambda_2_3>=0)&&(lambda_4>=0)){ //supersonic outlet{
		ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=		rho_v_y;
		ghost_state[3]= rho_energy;
			return 1;
	}else if((lambda_1<0)&&(lambda_2_3<0)&&(lambda_4<0)){//supersonic inlet
		ghost_state[0] =  rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=		rho_v_y_ext;
		ghost_state[3]= rho_energy_ext;
			return 2;			
	}else return 5;
}

void RiemannInvariants::get_du_du(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext,int bdry, int entry_j, double* dudu){


	double c = QuantityCalculator::calc_sound_speed(rho, rho_v_x, rho_v_y, rho_energy,kappa);
	double c_ext = QuantityCalculator::calc_sound_speed(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext,kappa);

if(bdry ==0){ //solid wall
				if(entry_j==0){
						dudu[0] = 1.;
						dudu[1] = 0.;
						dudu[2] = 0.;
						dudu[3] = 0.;
				}else if(entry_j==1){
						dudu[0] = 0.;
						dudu[1] = (n_y*n_y-n_x*n_x);
						dudu[2] = (-2*n_x*n_y);
						dudu[3] = 0.;
				}else if(entry_j==2){
						dudu[0] = 0.;
						dudu[1] = (-2*n_x*n_y);
						dudu[2] = (n_x*n_x-n_y*n_y);
						dudu[3] = 0.;
				}else if(entry_j==3){
						dudu[0] = 0.;
						dudu[1] = 0.;
						dudu[2] = 0.;
						dudu[3] = 1.;
				}else  throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no entry_j!");


	}else   throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no subsonic stream!");



}


void RiemannInvariants::get_free_stream(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, double t_x, double t_y,
					double* new_variables, 
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, int& boundary, bool solid ){

	double lambda_1 = RiemannInvariants::get_ev1(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_2_3 = RiemannInvariants::get_ev2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_4 = RiemannInvariants::get_ev4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);


	double rho_new= 0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
if(solid==true){  //solid wall
			rho_new = rho;
			rho_energy_new = rho_energy;
			rho_v_x_new = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
			rho_v_y_new = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
			boundary = 0;
	}else{
			if((lambda_1>=0)&&(lambda_2_3>=0)&&(lambda_4>=0)){ //supersonic outlet{

					rho_new= rho; 
					rho_v_x_new=rho_v_x; 	
					rho_v_y_new = rho_v_y;
					rho_energy_new = rho_energy;
					boundary = 1;
		}else if((lambda_1<0)&&(lambda_2_3<0)&&(lambda_4<0)){//supersonic inlet

					rho_new = rho_ext; 
					rho_v_x_new = rho_v_x_ext; 	
					rho_v_y_new = rho_v_y_ext;
					rho_energy_new = rho_energy_ext;
					boundary = 2;			
		}
	}

	new_variables[0]=  rho_new; 
	new_variables[1] = rho_v_x_new; 
	new_variables[2] = rho_v_y_new;
	new_variables[3] = rho_energy_new;
}


