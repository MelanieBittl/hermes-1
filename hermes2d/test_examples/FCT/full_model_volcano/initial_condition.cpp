#include "initial_condition.h"


//------------------- Initial condition ----------------
 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{ 
double pressure = pressure_init;
double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;
double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;  //barometrische hoehenformel

double rho = pressure/(temp_actual*R_gas);

		if(particle)	rho = rho_p*alpha_2;
		else rho*=(1.-alpha_2);	


return rho;

};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh, gamma,this->particle);
 };
//-------rho_v_x////Partikel geschwindigkeit != Gas geschwindigkeit!!!! sonst probleme bei Berechnung von source
 void CustomInitialCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.;
		dy= 0.;		
};

 double CustomInitialCondition_rho_v_x::value(double x, double y) const 
	{     


double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);

double v_x = 0;

double mach =mach_init;
if(!particle) v_x =	std::sqrt(gamma*pressure/rho)*mach;
//if(!particle) v_x*=vel_coeff;

if(particle) rho =alpha_2*rho_p;
else rho*=(1-alpha_2);	

return (rho*v_x);	
};

 Ord CustomInitialCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_x::clone() const {     
			return new CustomInitialCondition_rho_v_x(this->mesh, gamma,this->particle);
 };
//-------rho_v_y
 void CustomInitialCondition_rho_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.;
		dy= 0.;		
};

 double CustomInitialCondition_rho_v_y::value(double x, double y) const 
	{     

double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);
double v_y =0.; 

if(particle) rho =alpha_2*rho_p;
else rho*=(1.-alpha_2);		
	
	
return (rho*v_y);	
};

 Ord CustomInitialCondition_rho_v_y::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_y::clone() const {     
			return new CustomInitialCondition_rho_v_y(this->mesh, gamma,this->particle);
 };


//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;

};

 double CustomInitialCondition_e::value(double x, double y) const 
{  



double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;


double rho = pressure/(temp_actual*R_gas);

double v_x = 0;
double mach =mach_init;
 if(!particle) v_x =	std::sqrt(gamma*pressure/rho)*mach;
//if(!particle) v_x*=vel_coeff;

double v_y =0.; 

if(particle){
	rho =rho_p*alpha_2;
}else{
	rho*=(1.-alpha_2);	
	pressure *= (1.-alpha_2);
}
double rho_v_x = rho*v_x;
double rho_v_y = rho*v_y;

double temperature = temp_actual;
double energy = 0;
if(!particle) energy =QuantityCalculator::calc_energy(rho, rho_v_x ,rho_v_y, pressure, this->gamma);
else energy=(c_vp*temperature+0.5*(v_x*v_x+v_y*v_y))*rho;

	return energy;


};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,gamma,this->particle);

    }




//-----------------------------
//------bdry------------

//-------------------------------------



 void BoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho::value(double x, double y) const 
	{    			
double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);

	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
		if(particle)rho=rho_p*alpha;	
		else rho =(1.-alpha)*rho_g;
	}else{ 
		if(particle)	rho = rho_p*alpha_2;
		else rho*=(1-alpha_2);	
	
	}
return rho;
	
};

 Ord BoundaryCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho::clone() const {     
			return new BoundaryCondition_rho(this->mesh, gamma,this->particle);
 };

// v_x
 void BoundaryCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx = 0.; dy= 0.;


};

 double BoundaryCondition_rho_v_x::value(double x, double y) const 
	{   
double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) rho = rho_g;
double v_x =0.; 

double mach =mach_left;

if(!particle){
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) v_x=0;
else v_x =	std::sqrt(gamma*pressure/rho)*mach;
}

	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
		if(particle)rho=rho_p*alpha;	
		else rho =(1.-alpha)*rho_g;
	}else{ 
		if(particle)	rho = rho_p*alpha_2;
		else rho*=(1.-alpha_2);	
	
	}

return (rho*v_x);	


};

 Ord BoundaryCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_x::clone() const {     
			return new BoundaryCondition_rho_v_x(this->mesh,gamma,this->particle);
 };
//v_y
 void BoundaryCondition_rho_v_y::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx = 0.; dy= 0.;


};

 double BoundaryCondition_rho_v_y::value(double x, double y) const 
	{   
double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) rho = rho_g;
	//if(particle) rho = rho_p;

double v_y =0.; 
//double pressure = pressure_init;
//double mach =mach_volcano;
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
	/*		double factor = -1./22.5*(x-1.5)*(x-1.5)+1.;	
	pressure = pressure_coeff*pressure_init;
	v_y =std::sqrt(gamma*pressure/rho)*mach;
	if(!particle) v_y*=vel_coeff;
		v_y*=factor;*/
if(particle) v_y = v_y_p;
else v_y = v_y_g; 
}

	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
		if(particle)rho=rho_p*alpha;	
		else rho =(1.-alpha)*rho_g;
	}else{ 
		if(particle)	rho = rho_p*alpha_2;
		else rho*=(1.-alpha_2);	
	
	}

return (rho*v_y);	

};

 Ord BoundaryCondition_rho_v_y::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_y::clone() const {     
			return new BoundaryCondition_rho_v_y(this->mesh,gamma,this->particle);
 };
//bry_e
 void BoundaryCondition_rho_e::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx =0.;
		dy = 0;	
};

 double BoundaryCondition_rho_e::value(double x, double y) const 
	{  

double pressure = pressure_init;

double temp_actual = (288.15-0.0065*(y+3)*L)/T_infty;

double var = 1.-(0.0065*(y+3)*L/288.15);
double p_multiplied = std::pow(var,5.255);
pressure*=p_multiplied;

double rho = pressure/(temp_actual*R_gas);
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) rho = rho_g;
	//if(particle) rho = rho_p;



if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) pressure = pressure_coeff*pressure_init;

double v_x =0.; 
double mach =mach_left;
if(!particle){
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) v_x=0;
else v_x =	std::sqrt(gamma*pressure/rho)*mach;
}

double v_y =0.; 
//mach =mach_volcano;
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
	/*		double factor = -1./22.5*(x-1.5)*(x-1.5)+1.;	
	v_y =std::sqrt(gamma*pressure/rho)*mach;
	if(!particle) v_y*=vel_coeff;
 v_y*=factor;*/
	if(particle) v_y = v_y_p;
else v_y = v_y_g;

}

	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
		if(particle)rho=rho_p*alpha;	
		else{ rho =(1.-alpha)*rho_g;
				pressure *=(1.-alpha);
		}
	}else{ 
		if(particle)	rho = rho_p*alpha_2;
		else{ rho*=(1.-alpha_2);	
		pressure *=(1.-alpha_2);
		}
	
	}

double rho_v_x = rho*v_x;
double rho_v_y = rho*v_y;

double temperature = temp_actual;
if((y<=y_max)&&(x>=1.2)&&(x<=1.8)) temperature = temp;
double energy = 0;
if(!particle) energy =QuantityCalculator::calc_energy(rho, rho_v_x ,rho_v_y, pressure, this->gamma);
else energy=(c_vp*temperature+0.5*(v_x*v_x+v_y*v_y))*rho;

	return energy;
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,gamma, this->particle);
 };

