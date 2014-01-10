#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;


};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{    			

		return 1.;

};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh, kappa);
 };
//-------rho_v_x
 void CustomInitialCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx =0.;
		dy = 0.0;
		
};

 double CustomInitialCondition_rho_v_x::value(double x, double y) const 
	{     

return 1.;
			
};

 Ord CustomInitialCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_x::clone() const {     
			return new CustomInitialCondition_rho_v_x(this->mesh, kappa);
 };


//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;

};

 double CustomInitialCondition_e::value(double x, double y) const 
{  

	return QuantityCalculator::calc_energy(1.,1.,1.,1., this->kappa);


};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }

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
			return new BoundaryCondition_rho(this->mesh,this->time, kappa);
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
			return new BoundaryCondition_rho_v_x(this->mesh,this->time, kappa);
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

	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->kappa);
};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh,this->time, kappa);
 };

