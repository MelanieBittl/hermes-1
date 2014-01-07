#include "initial_condition.h"

//------------------- Initial condition ----------------
//----------density
 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
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
//----------rho_v_x_


 void CustomInitialCondition_v_x_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_v_x_rho::value(double x, double y) const 
	{       
			return 0.2;
			//return 1.;
};

 Ord CustomInitialCondition_v_x_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_v_x_rho::clone() const {     
			return new CustomInitialCondition_v_x_rho(this->mesh, kappa);
 };





//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;	
};

 double CustomInitialCondition_e::value(double x, double y) const 
{   


double pressure = 1./kappa;
double rho = 1.;
double rho_v_x = 0.2;
	return QuantityCalculator::calc_energy(rho, rho_v_x ,0.0, pressure, this->kappa);


};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }
