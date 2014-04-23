#include "euler_util.h"
#include "limits.h"
#include <limits>

// Calculates rho_energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double gamma)
{
  double to_return = pressure/(gamma - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
  return to_return;
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double rho_energy, double gamma)
{
  double to_return = (gamma - 1.0) * (rho_energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
  return to_return;
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double gamma)
{
  double to_return = std::sqrt(gamma * calc_pressure(rho, rho_v_x, rho_v_y, energy, gamma) / rho);
  return to_return;
}

  // Calculates enthalpy.
double QuantityCalculator::enthalpy(double rho, double rho_v_x, double rho_v_y, double rho_energy, double gamma){

double to_return= (rho_energy+ QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y,rho_energy,gamma))/rho;
return to_return;
}

//----------------------------------------------------------------------------------
//----------------------Filters.--------------------------------------
//----------------------------------------------------------------------------------
void MachNumberFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(gamma * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], gamma) / values.at(0)[i]);
}


void PressureFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (gamma - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}

void VelocityFilter_x::filter_fn(int n, Hermes::vector<double*> values, double* result)
{ 
  for (int i = 0; i < n; i++)
    result[i] = values.at(1)[i]/(values.at(0)[i]);
}
void VelocityFilter_y::filter_fn(int n, Hermes::vector<double*> values, double* result)
{ 
  for (int i = 0; i < n; i++)
    result[i] = values.at(2)[i]/(values.at(0)[i]);
}

void RadiusVelocityFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt( values.at(1)[i]*values.at(1)[i]+values.at(2)[i]*values.at(2)[i])/values.at(0)[i];
}

void EntropyFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], gamma) / p_ext)
      / Hermes::pow((values.at(0)[i] / rho_ext), gamma));
}

void AlphaFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{ 
  for (int i = 0; i < n; i++)
    result[i] = values.at(0)[i]/rho_p;
}

//-----------------------------------------------------
//-------------------RiemannInvariants-----------------------
//-----------------------------------------------------
double RiemannInvariants::get_w1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - 2.*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, gamma)/(gamma-1.));
}

double RiemannInvariants::get_w2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return QuantityCalculator::calc_pressure(rho, rho_v_x,rho_v_y, rho_energy, gamma)/Hermes::pow(rho,gamma);
}

double RiemannInvariants::get_w3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double t_x, double t_y){
		return  (rho_v_x*t_x/rho + rho_v_y*t_y/rho);
}

double RiemannInvariants::get_w4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +2*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, gamma)/(gamma-1));
}

bool RiemannInvariants::solid_wall(double rho, double rho_v_x, double rho_v_y, double n_x, double n_y){
	if(fabs(rho_v_x*n_x/rho + rho_v_y*n_y/rho)<1E-12) 
				return true;
	else false;

}


double RiemannInvariants::get_ev1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, gamma));
}

double RiemannInvariants::get_ev2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, gamma));
}

double RiemannInvariants::get_pressure(double w_1, double w_2, double w_3, double w_4){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return (c*c*rho/gamma);
}

double RiemannInvariants::get_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x){
	return (w_1+w_4)*n_x/2.+w_3*t_x;
}
double RiemannInvariants::get_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y){
	return (w_1+w_4)*n_y/2.+w_3*t_y;

}
double RiemannInvariants::get_speed_sound(double w_1, double w_4){
									return (gamma-1.)*(w_4-w_1)/4.;
}

double RiemannInvariants::get_rho(double w_1, double w_2, double w_3, double w_4){
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return  Hermes::pow(c*c/(gamma*w_2),(1/(gamma-1)));
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
	return (p/(gamma-1.)+rho*v_norm/2.);

}


void RiemannInvariants::get_ghost_state(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state)
{

		double w_1,w_2,w_3,w_4;	
	if( bdry==0){
		ghost_state[0] = rho;
		ghost_state[3] = rho_energy;
		ghost_state[1] = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
		ghost_state[2] = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
	}else if(bdry==1){ //supersonic outlet
		ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=		rho_v_y;
		ghost_state[3]= rho_energy;
	}else if(bdry==2){//supersonic inlet
		ghost_state[0] =  rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=		rho_v_y_ext;
		ghost_state[3]= rho_energy_ext;	
	}else if(bdry==3){//subsonic inlet
		w_2 = RiemannInvariants::get_w2(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, n_x, n_y);
		w_3 = RiemannInvariants::get_w3(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, t_x, t_y);
		w_4 = RiemannInvariants::get_w4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
 		double c = QuantityCalculator::calc_sound_speed(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext,gamma);
		w_1 = w_4 - 4.*c/(gamma-1);

		ghost_state[0]= get_rho(w_1, w_2, w_3, w_4);
		ghost_state[1] =get_rho_v_x(w_1, w_2, w_3, w_4, n_x,  t_x);
		ghost_state[2] =get_rho_v_y(w_1, w_2, w_3, w_4, n_y,  t_y);
		ghost_state[3] =get_energy(w_1, w_2, w_3, w_4,  n_x,  n_y,  t_x,  t_y);
		
	}else if(bdry==4){//subsonic outlet
		w_2 = RiemannInvariants::get_w2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		w_3 = RiemannInvariants::get_w3(rho, rho_v_x, rho_v_y, rho_energy, t_x, t_y);
		w_4 = RiemannInvariants::get_w4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		double pressure_out =QuantityCalculator::calc_pressure(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, gamma);
		double pressure_in = QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
		w_1 = w_4 - 4./(gamma-1)*std::sqrt(gamma*pressure_out/rho*std::pow(pressure_in/pressure_out, 1./gamma));

		ghost_state[0]= get_rho(w_1, w_2, w_3, w_4);
		ghost_state[1] =get_rho_v_x(w_1, w_2, w_3, w_4, n_x,  t_x);
		ghost_state[2] =get_rho_v_y(w_1, w_2, w_3, w_4, n_y,  t_y);
		ghost_state[3] =get_energy(w_1, w_2, w_3, w_4,  n_x,  n_y,  t_x,  t_y);

	}
}

void RiemannInvariants::get_ghost_state_p(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state)
{
	if((bdry ==3)||(bdry ==2)) //inlet
	{	ghost_state[0] = rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=	 rho_v_y_ext;
		ghost_state[3]=  rho_energy_ext;	
	}else if((bdry ==1)||(bdry ==4)) //outlet
	{	ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=	 rho_v_y;
		ghost_state[3]=  rho_energy;	
	}

}



void RiemannInvariants::get_du_du(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext,int bdry, int entry_j, double* dudu){

double w_1,w_2,w_3,w_4;	
	double c = QuantityCalculator::calc_sound_speed(rho, rho_v_x, rho_v_y, rho_energy,gamma);
	double c_ext = QuantityCalculator::calc_sound_speed(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext,gamma);

	if(bdry ==3){//subsonic inlet
		w_2 = RiemannInvariants::get_w2(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, n_x, n_y);
		w_3 = RiemannInvariants::get_w3(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, t_x, t_y);
		w_4 = RiemannInvariants::get_w4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);		
 		w_1 = w_4 - 4.*c_ext/(gamma-1);
		double rho_ghost = get_rho(w_1, w_2, w_3, w_4);
		double v_x_ghost = get_v_x(w_1, w_2, w_3, w_4, n_x,  t_x);
		double v_y_ghost = get_v_y(w_1, w_2, w_3, w_4, n_y,  t_y);

			double dw4_du;
			if(entry_j==0){
					dw4_du= -1./rho*(rho_v_x*n_x+rho_v_y*n_y)/rho+ gamma/c*(-rho_energy/(rho*rho)+(rho_v_x*rho_v_x+rho_v_y*rho_v_y)/(rho*rho*rho));		
			}else if(entry_j==1){
					dw4_du = n_x/rho-gamma*rho_v_x/(rho*rho*c);
			}else if(entry_j==2){
					dw4_du = n_y/rho-gamma*rho_v_y/(rho*rho*c);
			}else if(entry_j==3){
					dw4_du = gamma/(c*rho);
			}else  throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no entry_j!");
			double drho_du = rho_ghost/(2.*c_ext)*dw4_du;

				dudu[0] = drho_du;
				dudu[1] = (v_x_ghost* drho_du + 0.5* n_x* rho_ghost*dw4_du);
				dudu[2] = (v_y_ghost * drho_du + 0.5* n_y* rho_ghost*dw4_du);
				double dp_du = ((gamma-1.)*0.5*rho_ghost*c_ext*dw4_du + c_ext*c_ext*drho_du)/gamma;
				dudu[3] = (1./(gamma-1.)*dp_du+ (v_x_ghost*v_x_ghost+v_y_ghost*v_y_ghost)*0.5*drho_du+0.25*rho_ghost*(w_1+w_4)*dw4_du);			

	}else if(bdry ==4){//outlet
		double du_du,dv_du,dp_du, dw4_du, dw2_du, dw3_du;
		w_2 = RiemannInvariants::get_w2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		w_3 = RiemannInvariants::get_w3(rho, rho_v_x, rho_v_y, rho_energy, t_x, t_y);
		w_4 = RiemannInvariants::get_w4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		double pressure_out =QuantityCalculator::calc_pressure(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, gamma);
		double pressure_in = QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
		w_1 = w_4 - 4./(gamma-1)*std::sqrt(gamma*pressure_out/rho*std::pow(pressure_in/pressure_out, 1./gamma));
		double rho_ghost = std::pow(c_ext*c_ext/(gamma*w_2), 1./(gamma-1.));
		double v_x_ghost = (0.5*n_x*(w_1+w_4)+ t_x*w_3);
		double v_y_ghost = (0.5*n_y*(w_1+w_4)+ t_y*w_3);

			if(entry_j==0){
					dw4_du = -1./rho*(rho_v_x*n_x+rho_v_y*n_y)/rho+ gamma/c*(-rho_energy/(rho*rho)+(rho_v_x*rho_v_x+rho_v_y*rho_v_y)/(rho*rho*rho));	
					dw2_du = (gamma-1.)*0.5/std::pow(rho,gamma)*(rho_v_x*rho_v_x+rho_v_y*rho_v_y)/(rho*rho)- gamma/rho*	w_2;
					dw3_du = 1./rho*(rho_v_x/rho*n_y - rho_v_y/rho*n_x);
			}else if(entry_j==1){
					dw4_du = n_x/rho-gamma*rho_v_x/(rho*rho*c);
					dw2_du = (1-gamma)/std::pow(rho,gamma)*rho_v_x/rho;
					dw3_du = -n_y/rho;
			}else if(entry_j==2){
					dw4_du = n_y/rho-gamma*rho_v_y/(rho*rho*c);
					dw2_du = (1.-gamma)/std::pow(rho,gamma)*rho_v_y/rho;
					dw3_du = n_x/rho;
			}else if(entry_j==3){
					dw4_du = gamma/(c*rho);
					dw2_du = (gamma-1.)/std::pow(rho,gamma);
					dw3_du = 0.;
			}else  throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no entry_j!");
			double drho_du = rho_ghost/(2.*c_ext)*dw4_du+ rho_ghost/((1.-gamma)*w_2)*dw2_du;
				dudu[0] =  drho_du;
				du_du= 0.5*n_x*dw4_du-n_y*dw3_du;
				dv_du = 0.5*n_y*dw4_du+n_x*dw3_du;
 				dudu[1] =  (rho_ghost*du_du + v_x_ghost*drho_du);
 				dudu[2] =  (rho_ghost*dv_du + v_y_ghost*drho_du);

				dv_du= 2.*w_3*dw3_du+0.5*(w_1+w_4)*dw4_du;
				dp_du = (gamma-1.)*0.5/gamma*rho_ghost*c_ext*dw4_du+1./gamma*c_ext*c_ext*drho_du;		
				dudu[3] =  (1./(1.-gamma)*dp_du+ 0.5*rho_ghost*dv_du + (v_x_ghost*v_x_ghost+v_y_ghost*v_y_ghost)/(2.)*drho_du );			
	
	}else   throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no subsonic stream bdry = %i!", bdry);



}

