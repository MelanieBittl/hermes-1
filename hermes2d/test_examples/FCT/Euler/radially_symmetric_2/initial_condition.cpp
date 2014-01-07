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

 double CustomInitialCondition_rho::value(double x, double y) const {       

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
			return new CustomInitialCondition_rho(this->mesh);
 };

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
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

 double CustomInitialCondition_e::value(double x, double y) const {       
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

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->kappa);

};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(2);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }







 void CustomInitialCondition_p::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_p::value(double x, double y) const {       
			if((x*x+y*y)>=(0.13*0.13)) 0.1;
			else			15.0;
};

 Ord CustomInitialCondition_p::ord(double x, double y)   const {
      return Ord(2);
};

 MeshFunction<double>* CustomInitialCondition_p::clone() const
    {
return new CustomInitialCondition_p(this->mesh);

    }
