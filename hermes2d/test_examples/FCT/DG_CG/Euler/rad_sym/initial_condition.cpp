#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const {       
			if((x*x+y*y)>=(0.13*0.13)) return 1.0;
			else			return 2.;

};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(2);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh);
 };

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			if((x*x+y*y)>=(0.13*0.13)) return QuantityCalculator::calc_energy(1.0, 0.0 ,0.0, 1.0, kappa);
			else			return QuantityCalculator::calc_energy(2.0, 0.0 ,0.0, 15.0, kappa);
};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(2);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }
