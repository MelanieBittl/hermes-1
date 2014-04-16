#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{ 
		
double rho = 1.; 

if(particle){rho = 0.1;
}
return rho;

};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh, gamma,this->particle);
 };
//-------rho_v_x
 void CustomInitialCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.;
		dy= 0.;		
};

 double CustomInitialCondition_rho_v_x::value(double x, double y) const 
	{     
double rho = 1.; 

double pressure = std::pow(10,1);
double mach = 0.01;

double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){
 rho = 0.1; 
	v_x = 0;
}

return (rho*v_x);	
};

 Ord CustomInitialCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_x::clone() const {     
			return new CustomInitialCondition_rho_v_x(this->mesh, gamma,this->particle);
 };


//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;

};

 double CustomInitialCondition_e::value(double x, double y) const 
{  

double rho = 1.; 

double pressure = std::pow(10,1);
double mach = 0.01;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){ v_x = 0; pressure = 0;
 rho = 0.1; 
}

double rho_v_x = rho*v_x;

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->gamma);

};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,gamma,this->particle);

    }




//-----------------------------
//------bdry
 void BoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho::value(double x, double y) const 
	{    			
		
double rho = 1.;
if(particle){rho = 10.; if(x>1) rho = 0.1;}
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
double rho = 1.;

double pressure = std::pow(10,1);

double mach = 0.25-y*y;
if(x>1) mach = 0.01;

double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){rho = 10.; if(x>1){ rho = 0.1;v_x = 0;}}

return (rho*v_x);	

};

 Ord BoundaryCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_x::clone() const {     
			return new BoundaryCondition_rho_v_x(this->mesh,gamma,this->particle);
 };


//bry_e
 void BoundaryCondition_rho_e::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx =0.;
		dy = 0;	

};

 double BoundaryCondition_rho_e::value(double x, double y) const 
	{    			

double rho = 1.;

double pressure = std::pow(10,1);
double mach = -y*y +0.25;
if(x>1) mach = 0.01;

double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){rho = 10.; if(x>1){ rho = 0.1;v_x = 0;pressure = 0;}}
double rho_v_x = rho*v_x;

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->gamma);
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,gamma, this->particle);
 };

