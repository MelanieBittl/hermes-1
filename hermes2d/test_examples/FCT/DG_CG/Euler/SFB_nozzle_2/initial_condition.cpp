#include "initial_condition.h"
const double y_1 = 10.3875; //13.45;
const double y_2 = 17.3125; //14.25;
//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{    			


double rho = 1.1881;

if((y >y_1)&&(y<y_2)) rho =1.778;

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


double v_x =0.;



double rho = 1.1881;
if((y >y_1)&&(y<y_2)){ rho =1.778;v_x = 442.12;}
double pressure = 100000.;



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


double v_y =0.;
double rho = 1.1881;
if((y >y_1)&&(y<y_2)) rho =1.778;
double pressure = 100000.;
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
		dx =0.;	dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const 
{  


double v_x= 0.; 

double rho = 1.1881;
if((y >y_1)&&(y<y_2)){ rho =1.778;v_x = 442.12;}
double pressure = 100000.;
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


double rho = 1.778;
double pressure = 100000.;




return rho;

};

 Ord BoundaryCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho::clone() const {     
			return new BoundaryCondition_rho(this->mesh, gamma,this->particle);
 };

//------------ v_x
 void BoundaryCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx = 0.; dy= 0.;


};

 double BoundaryCondition_rho_v_x::value(double x, double y) const 
	{    
		

double rho = 1.778;
double pressure = 100000.;

double v_x = 442.12;

return (rho*v_x);	

};

 Ord BoundaryCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_x::clone() const {     
			return new BoundaryCondition_rho_v_x(this->mesh,gamma,this->particle);
 };
 //------vy
 void BoundaryCondition_rho_v_y::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx = 0.; dy= 0.;


};

 double BoundaryCondition_rho_v_y::value(double x, double y) const 
	{    
		
double rho = 1.778;
double pressure = 100000.;

double v_x = 442.12;

return 0;

};

 Ord BoundaryCondition_rho_v_y::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_y::clone() const {     
			return new BoundaryCondition_rho_v_y(this->mesh,gamma,this->particle);
 };



//-------------bry_e
 void BoundaryCondition_rho_e::derivatives(double x, double y, double& dx, double& dy) const 
{ 
		dx =0.;
		dy = 0;	

};

 double BoundaryCondition_rho_e::value(double x, double y) const 
	{    			

double rho = 1.778;
double pressure = 100000.;

double v_x = 442.12;

double rho_v_x = rho*v_x;
double rho_v_y = rho*0.;


	return QuantityCalculator::calc_energy(rho, rho_v_x ,rho_v_y, pressure, this->gamma);


	
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,gamma, this->particle);
 };

