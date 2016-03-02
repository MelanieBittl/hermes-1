#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const {       
return 1.;

};

 Ord CustomInitialCondition_rho::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh);
 };
 
 
 
  void CustomInitialCondition_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_v_x::value(double x, double y) const {       
return 0.67;

};

 Ord CustomInitialCondition_v_x::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition_v_x::clone() const {     
			return new CustomInitialCondition_v_x(this->mesh);
 };

 void CustomInitialCondition_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_v_y::value(double x, double y) const {       
return 0;

};

 Ord CustomInitialCondition_v_y::ord(double x, double y)  const {
      return Ord(1);
};
 
 MeshFunction<double>* CustomInitialCondition_v_y::clone() const {     
			return new CustomInitialCondition_v_y(this->mesh);
 }; 
 

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			return QuantityCalculator::calc_energy(1., 0.67 ,0.0, 1.0/kappa, kappa);

};

 Ord CustomInitialCondition_e::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }

