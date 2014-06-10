#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{ 
		
double rho = 1.; 

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
double rho = 1.; 
double pressure = 1.;
double mach = 1.;
if((x<=2)&&(x>=1.5)&&(y<-0.5)) return 0;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
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
double rho = 1.; 
double pressure = 1.;
double mach = 1.;
if((x<=2)&&(x>=1.5)&&(y<-0.5)) mach = 0.8;
else mach =  0.;
double v_y = std::sqrt(gamma*pressure/rho)*mach;
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

double rho = 1.; 
double pressure = 1.; 
double mach =1.;
if((x<=2)&&(x>=1.5)&&(y<-0.5)) mach = 0.8;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){ pressure = 0; rho=0.1; v_x *=0.9;
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
//------bdry------------

//-------------------------------------



 void BoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho::value(double x, double y) const 
	{    			
	/*	
double rho = 6.; 

if(particle){ if(x>1) rho = 1.9-x*0.9; else rho=1;}
return rho;*/
	double rho = 1.; 

if(particle){ rho=0.1;
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
		
		double rho = 1.; 
double pressure = 10.;

double mach = 1.;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){ rho=0.1; v_x *=0.9;
}
if(y<-1) return 0;
return (rho*v_x);	
/*
double rho = 6.; 

double pressure = 1;
double mach = 0;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){rho =1;v_x +=0.1;}

return (rho*v_x);	*/

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
		
		double rho = 1.; 
double pressure = 10.;
if(y>-2) return 0;
double mach = 1.;
double v_y = std::sqrt(gamma*pressure/rho)*mach;
if(particle){ rho=0.1; v_y *=0.9;
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

/*double rho = 6.; 

double pressure = 1;
double mach = 0;
double v_x = std::sqrt(gamma*pressure/rho)*mach;

if(particle){rho = 1;pressure = 0;v_x +=0.1;}
double rho_v_x = rho*v_x;*/

double rho = 1.; 
double pressure = 10.; 
double mach =1;
if(y<-1) return 100;
double v_x = std::sqrt(gamma*pressure/rho)*mach;
if(particle){ pressure = 0; rho=0.1; v_x *=0.9;
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

