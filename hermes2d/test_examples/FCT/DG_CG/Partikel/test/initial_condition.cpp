#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
		
  
	double radius = 0.;
     //hump
	double x_0 =0.;
	double y_0= 0.;	
	radius = (2.) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) {		
		dx = -std::sin(radius*PI)/4.0*(PI/(0.5 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(x-x_0);
		dy = -std::sin(radius*PI)/4.0*(PI/(0.5 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(y-y_0);	
	}

};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{    			

		//return 1.;

    double result = 0.0;
 	double radius;
       //hump
	double x_0 =0.;
	double y_0= 0.;	
	radius = 2.* std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) { 
		 result = (1.0+ std::cos(PI*radius))/4.0; 			
		
	}
return result+1.;	
};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh, gamma);
 };
//-------rho_v_x
 void CustomInitialCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx =0.;
		dy = 0.0;
		
};

 double CustomInitialCondition_rho_v_x::value(double x, double y) const 
	{     

//	return 3.5;
return 0.;
			
};

 Ord CustomInitialCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_x::clone() const {     
			return new CustomInitialCondition_rho_v_x(this->mesh, gamma);
 };


//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
		
		dx =0.;
		dy = 0.0;
		
  
	double radius = 0.; double result = 0.;
     //hump
	double x_0 =0.;
	double y_0= 0.;	
	radius = (2.) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) {		
		dx = -std::sin(radius*PI)/4.0*(PI/(0.5 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(x-x_0);
		dy = -std::sin(radius*PI)/4.0*(PI/(0.5 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(y-y_0);
result = (1.0+ std::cos(PI*radius))/4.0; 	
	}

dx*=1./std::pow(result+1,2.0);
dy*=1./std::pow(result+1,2.0);
};

 double CustomInitialCondition_e::value(double x, double y) const 
{  

 	double radius; double result = 0.;
       //hump
	double x_0 =0.;
	double y_0= 0.;	

	radius = 2.*std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) { 
		 result = (1.0+ std::cos(PI*radius))/4.0; 
					
		
	}

double pressure = 1/(result+1);  
double  rho= result+1.;

double rho_v_x = 0.;

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->gamma);


};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,gamma);

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
//return 1.;
		return Hermes::sin(time)+1.;
};

 Ord BoundaryCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho::clone() const {     
			return new BoundaryCondition_rho(this->mesh,this->time, gamma);
 };

// v_x
 void BoundaryCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho_v_x::value(double x, double y) const 
	{    			
		//return 3.5;
		return (Hermes::sin(time)+1.)*3.5;
};

 Ord BoundaryCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_x::clone() const {     
			return new BoundaryCondition_rho_v_x(this->mesh,this->time, gamma);
 };


//bry_e
 void BoundaryCondition_rho_e::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho_e::value(double x, double y) const 
	{    			

double pressure =1.;// (Hermes::sin(time)+1.);
double  rho=  (Hermes::sin(time)+1.);

double rho_v_x = rho*3.5;

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->gamma);
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,this->time, gamma);
 };

