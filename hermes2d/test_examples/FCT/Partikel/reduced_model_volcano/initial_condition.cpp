#include "initial_condition.h"

double y_max =1.;

double mach_volcano = 1.;
double mach_left = 0.01;

double pressure_init = 101325; //std::pow(10,5); //1MPa= 10^6, 0.101325 MPa = atmospheric pressure
double pressure_coeff= 2.;
double temp = 1200;
double temp_init = 298;
double R_gas=287.;
double rho_p = 2300.; 
double rho_g = pressure_init*pressure_coeff/(temp*R_gas); //p = rho*R*T
double alpha = 0.01; //particle ash




//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{ 
		
double rho = pressure_init/(temp_init*R_gas);
if(particle) rho =0;

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

double rho = pressure_init/(temp_init*R_gas);

double v_x =0.; 

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

double rho = pressure_init/(temp_init*R_gas);
double v_y =0;
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
double rho = pressure_init/(temp_init*R_gas);
double v_x = 0;
double v_y =0.; 
double pressure = pressure_init;
double rho_v_x = rho*v_x;
double rho_v_y = rho*v_y;


	return QuantityCalculator::calc_energy(rho, rho_v_x ,rho_v_y, pressure, this->gamma);


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
double rho = pressure_init/(temp_init*R_gas);
	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
		if(particle)rho=rho_p*alpha;	
		else rho =(1-alpha)*rho_g;
	}else{ 
		if(particle) rho = 0;	
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
double rho = pressure_init/(temp_init*R_gas);
	if((y<=y_max)&&(x>=1.2)&&(x<=1.8)){
			rho =(1-alpha)*rho_g;
	}
double v_x =0.; 
double pressure = pressure_init;

double mach =mach_left;

if((y<=1.)&&(x>=1.2)&&(x<=1.8)) v_x=0;
else v_x =	std::sqrt(gamma*pressure/rho)*mach;




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
double rho = pressure_init/(temp_init*R_gas);
	if((y<=y_max)&&(x>=1.2)&&(x<=1.8))
				rho =(1-alpha)*rho_g;

double mach =mach_left;
double v_x =0.; 

double pressure =pressure_init;
if((y<=1.)&&(x>=1.2)&&(x<=1.8)) 	pressure = pressure_coeff*pressure_init;

if((y<=1.)&&(x>=1.2)&&(x<=1.8)) v_x=0;
else v_x =	std::sqrt(gamma*pressure/rho)*mach;


double v_y =0.;
mach =mach_volcano;
if((y<=1.)&&(x>=1.2)&&(x<=1.8)){
	v_y=std::sqrt(gamma*pressure/rho)*mach;
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
		
double rho = pressure_init/(temp_init*R_gas);
	if((y<=y_max)&&(x>=1.2)&&(x<=1.8))
				rho =(1-alpha)*rho_g;
	
double v_x =0.; 
double pressure =pressure_init;
if((y<=1.)&&(x>=1.2)&&(x<=1.8)) 	pressure = pressure_coeff*pressure_init;

double v_y =0.;
double mach =mach_volcano;
if((y<=1.)&&(x>=1.2)&&(x<=1.8)){
	v_y=std::sqrt(gamma*pressure/rho)*mach;
}

double rho_v_x = rho*v_x;
double rho_v_y = rho*v_y;


	return QuantityCalculator::calc_energy(rho, rho_v_x ,rho_v_y, pressure, this->gamma);
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,gamma, this->particle);
 };

