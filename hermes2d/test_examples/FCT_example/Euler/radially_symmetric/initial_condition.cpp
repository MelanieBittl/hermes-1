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

 Ord CustomInitialCondition_rho::ord(Ord x, Ord y) const {
      return Ord(2);
};

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			if((x*x+y*y)>=(0.13*0.13)) return QuantityCalculator::calc_energy(1.0, 0.0 ,0.0, 1.0, kappa);
			else			return QuantityCalculator::calc_energy(2.0, 0.0 ,0.0, 15.0, kappa);
};

 Ord CustomInitialCondition_e::ord(Ord x, Ord y) const {
      return Ord(2);
};


 void CustomInitialCondition_p::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_p::value(double x, double y) const {       
			if((x*x+y*y)>=(0.13*0.13)) 0.1;
			else			15.0;
};

 Ord CustomInitialCondition_p::ord(Ord x, Ord y) const {
      return Ord(2);
};


