#include "initial_condition.h"

//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_rho::value(double x, double y) const {       
			if(x< (1./6.+y/std::sqrt(3.))) return 8.0;
			else			return 1.4;

};

 Ord CustomInitialCondition_rho::ord(Ord x, Ord y) const {
      return Ord(2);
};

 void CustomInitialCondition_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_v_x::value(double x, double y) const {       
			if(x< (1./6.+y/std::sqrt(3.))) return 8.25*std::cos(M_PI/6.)*8.0;
			else			return 0.0;

};

 Ord CustomInitialCondition_v_x::ord(Ord x, Ord y) const {
      return Ord(2);
};
 void CustomInitialCondition_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_v_y::value(double x, double y) const {       
			if(x< (1./6.+y/std::sqrt(3.))) return -8.25*std::sin(M_PI/6.)*8.0;
			else			return 0.0;

};

 Ord CustomInitialCondition_v_y::ord(Ord x, Ord y) const {
      return Ord(2);
};



 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			if(x< (1./6.+y/std::sqrt(3.))) return QuantityCalculator::calc_energy(8.0, 8.25*std::cos(M_PI/6.)*8.0 ,-8.25*std::sin(M_PI/6.)*8.0, 116.5, kappa);
			else			return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);
};

 Ord CustomInitialCondition_e::ord(Ord x, Ord y) const {
      return Ord(2);
};



//Boundary-Condition


 void CustomBoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_rho::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return 8.0;
			else			return 1.4;



};

 Ord CustomBoundaryCondition_rho::ord(Ord x, Ord y) const {
      return Ord(2);
};

 void CustomBoundaryCondition_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_v_x::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return 8.25*std::cos(M_PI/6.)*8.0;
			else			return 0.0;

};

 Ord CustomBoundaryCondition_v_x::ord(Ord x, Ord y) const {
      return Ord(2);
};
 void CustomBoundaryCondition_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_v_y::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return -8.25*std::sin(M_PI/6.)*8.0;
			else			return 0.0;

};

 Ord CustomBoundaryCondition_v_y::ord(Ord x, Ord y) const {
      return Ord(2);
};



 void CustomBoundaryCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_e::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return QuantityCalculator::calc_energy(8.0, 8.25*std::cos(M_PI/6.)*8.0 ,-8.25*std::sin(M_PI/6.)*8.0, 116.5, kappa);
			else			return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);
};

 Ord CustomBoundaryCondition_e::ord(Ord x, Ord y) const {
      return Ord(2);
};

