#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =-1;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{    			

/*double  rho= 2.5;
if(x<0.5) rho = 3.0;
else if(x<1) rho = 3.0 -(x-0.5);*/
double rho = 1.;
if(x>1.) rho = 0.125;

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

/*double  rho=2.5;
if(x<0.5) rho = 3.0;
else if(x<1) rho = 3.0 -(x-0.5);


double pressure = 2.5;
if(x<0.5) pressure = 3.0;
else if (x<1) pressure = 3.0 -(x-0.5);*/

double v_x =0.; //std::sqrt(gamma*pressure/3.);

double rho = 1.;
if(x>1.) rho = 0.125;

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
		dx =0.;	dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const 
{  
/*
double  rho= 2.5;
if(x<0.5) rho = 3.0;
else if(x<1)  rho = 3.0 -(x-0.5);

double pressure = 2.5;
if(x<0.5) pressure = 3.0;
else if(x<1) pressure = 3.0 -(x-0.5);*/

double v_x= 0.; //=std::sqrt(gamma*pressure/3.);

double rho = 1.;
if(x>1.) rho = 0.125;

double pressure = 1.;
if(x>1.) pressure = 0.1;

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

/*double  rho= 3.0;

if(x>1) rho= 2.5;*/


double rho = 1.;
if(x>1.) rho = 0.125;

double pressure = 1.;
if(x>1.) pressure = 0.1;

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
		
/*
double  rho= 3.0;
double pressure = 3.;
double max = 1.5*std::sqrt(gamma*pressure/rho);
double v_x =max;//-Hermes::sqr(2.*y)*(max-0.5)+max;

if(x>1){ rho= 2.5; pressure= 2.5;
v_x = 0;//std::sqrt(gamma*pressure/3.);
}*/

double rho = 1.;
if(x>1.) rho = 0.125;

double pressure = 1.;
if(x>1.) pressure = 0.1;

double v_x = 0.;

if(x>1){ //rho= 2.5; pressure= 2.5;
v_x = 0;//std::sqrt(gamma*pressure/1.5);
}
double rho_v_x = rho*v_x;

return rho_v_x;

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

/*
double  rho= 3.0;
double pressure = 3.;
double max = 1.5*std::sqrt(gamma*pressure/rho);
double v_x =max;//-Hermes::sqr(2.*y)*(max-0.5)+max;

if(x>1){ rho= 2.5; pressure= 2.5;
v_x = 0;//std::sqrt(gamma*pressure/3.);
}*/



double rho = 1.;
if(x>1.) rho = 0.125;

double pressure = 1.;
if(x>1.) pressure = 0.1;

double v_x = 0.;

if(x>1){ //rho= 2.5; pressure= 2.5;
v_x = 0;//std::sqrt(gamma*pressure/1.5);
}
double rho_v_x = rho*v_x;


	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->gamma);


	
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,gamma, this->particle);
 };

