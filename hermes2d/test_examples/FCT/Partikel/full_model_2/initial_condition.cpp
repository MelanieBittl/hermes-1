#include "initial_condition.h"
const double y_1 = 12.35; 
const double y_2 = 15.35; 
const double multiplier = 20;

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{ 
		
double rho = 1.1881;

if(particle){ rho=0.1;
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
double v_x =0.;
double rho = 1.1881;
if(particle){ rho=0.1; v_x *=0.9;
}
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
if(particle){ rho=0.1; v_y *=0.9;
}
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

double v_x= 0.; 

double rho = 1.1881;double pressure = 100000.;
if(particle){ rho=0.1; v_x *=0.9; pressure =0.;}

double rho_v_x = rho*v_x;

double rho_v_y = 0;


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

double v_x= 0.; 
double rho = 1.1881;
if(particle){ rho=0.1;}
if((y >y_1)&&(y<y_2)){ rho =1.778;v_x = 442.12;}

if(x>16.7){ rho =1.778;}


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
		
double v_x= 0.; 
double rho = 1.1881;
double pressure = 100000.;
if((y >=y_1)&&(y<=y_2)){ rho =1.778;v_x = 442.12;}
if(x>16.7){ rho =1.778;v_x = 16.8*multiplier; }

if(particle){ rho=0.1; v_x *=0.9;
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
		
		
double v_y= 0.; 
double v_x=0.;
double rho = 1.1881;
double pressure = 100000.;
if((y >=y_1)&&(y<=y_2)){ rho =1.778;v_x = 442.12;}
if(x>16.7){ rho =1.778;v_y = 21.1*multiplier;
if(y>20) v_y*=-1;}

if(particle){ rho=0.1; v_y *=0.9;
}
return rho*v_y;


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

double v_x= 0.;  double v_y =0.;
double rho = 1.1881;
double pressure = 100000.;
if((y >=y_1)&&(y<=y_2)){ rho =1.778;v_x = 442.12;}
if(x>16.7){ rho =1.778;v_y = 21.1*multiplier;v_x = 16.8*multiplier;
if(y>20) v_y*=-1;};

if(particle){ pressure = 0; rho=0.1; v_x *=0.9;
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

