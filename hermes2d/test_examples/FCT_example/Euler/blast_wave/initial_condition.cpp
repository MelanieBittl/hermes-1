#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const {       
				return 1.0;

};

 Ord CustomInitialCondition_rho::ord(Ord x, Ord y) const {
      return Ord(2);
};

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			if(x>=0.90) return QuantityCalculator::calc_energy(1.0, 0.0 ,0.0, 100., kappa);
			else if(x<0.1) return QuantityCalculator::calc_energy(1.0, 0.0 ,0.0, 1000.0, kappa);
			else			return QuantityCalculator::calc_energy(1.0, 0.0 ,0.0, 0.01, kappa);
};

 Ord CustomInitialCondition_e::ord(Ord x, Ord y) const {
      return Ord(2);
};





