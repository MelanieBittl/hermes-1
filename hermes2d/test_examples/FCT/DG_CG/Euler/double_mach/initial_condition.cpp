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
			if(x< (1./6.+y/std::sqrt(3.))) return 8.25*std::sqrt(3)/2.*8.0;
			else			return 0.0;

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
			if(x< (1./6.+y/std::sqrt(3.))) return -8.25*0.5*8.0;
			else			return 0.0;

};

 Ord CustomInitialCondition_v_y::ord(double x, double y)  const {
      return Ord(10);
};
 
 MeshFunction<double>* CustomInitialCondition_v_y::clone() const {     
			return new CustomInitialCondition_v_y(this->mesh);
 }; 
 

 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomInitialCondition_e::value(double x, double y) const {       
			if(x< (1./6.+y/std::sqrt(3.))) return QuantityCalculator::calc_energy(8.0, 8.25*std::sqrt(3)/2.*8.0 ,-8.25*0.5*8.0, 116.5, kappa);
			else			return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);

};

 Ord CustomInitialCondition_e::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,kappa);

    }

//------------------- Boundary ----------------

 void CustomBoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_rho::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return 8.0;
			else			return 1.4;

};

 Ord CustomBoundaryCondition_rho::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomBoundaryCondition_rho::clone() const {     
			return new CustomBoundaryCondition_rho(this->mesh, this->time);
 };
 
 
 
  void CustomBoundaryCondition_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_v_x::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return 8.25*std::sqrt(3)/2.*8.0;
			else			return 0.0;

};

 Ord CustomBoundaryCondition_v_x::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomBoundaryCondition_v_x::clone() const {     
			return new CustomBoundaryCondition_v_x(this->mesh, this->time);
 };

 void CustomBoundaryCondition_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_v_y::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return -8.25*0.5*8.0;
			else			return 0.0;

};

 Ord CustomBoundaryCondition_v_y::ord(double x, double y)  const {
      return Ord(10);
};
 
 MeshFunction<double>* CustomBoundaryCondition_v_y::clone() const {     
			return new CustomBoundaryCondition_v_y(this->mesh, this->time);
 }; 
 

 void CustomBoundaryCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
		dx = 0.0;
		dy = 0.0;
};

 double CustomBoundaryCondition_e::value(double x, double y) const {       
			if(x< (1./6.+(1+20*time)/std::sqrt(3.))) return QuantityCalculator::calc_energy(8.0, 8.25*std::sqrt(3)/2.*8.0 ,-8.25*0.5*8.0, 116.5, kappa);
			else			return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);

};

 Ord CustomBoundaryCondition_e::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomBoundaryCondition_e::clone() const
    {
return new CustomBoundaryCondition_e(this->mesh,kappa, this->time);

    }

