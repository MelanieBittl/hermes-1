#include "initial_condition.h"

double calculate_a(double x, double y, double kappa)
{
	double q = 0.5;	
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
double rho, q_2, J, f_a, dJ,drho,dq_2,df_a, a_old;
	for(int i = 0; i<100;i++)
	{
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
			J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));	
		 f_a = std::pow(x+J*0.5, 2) + std::pow(y,2) - 1./(4*rho*rho*q_2*q_2);
		if(std::fabs(f_a)<1e-10) break;
		 dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+std::log(std::exp(1))/(1-a*a));
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 df_a= 	(x+0.5*J)*dJ + 0.25*(1./std::pow(q_2,2)*2*drho/std::pow(rho,3)+2*dq_2/(rho*rho*std::pow(q_2,3)));

		 a_old = a;
		a = a_old - f_a/df_a;

	}

return a; 
}


//------------------- Initial condition ----------------

 void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;

};

 double CustomInitialCondition_rho::value(double x, double y) const 
	{    			
	double a = calculate_a(x,y, this->kappa)	; 
	double	rho = std::pow(a, (2./(this->kappa-1.)));

return rho;	
};

 Ord CustomInitialCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho::clone() const {     
			return new CustomInitialCondition_rho(this->mesh, this->kappa);
 };
//-------rho_v_x
 void CustomInitialCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
		dx =0.;
		dy = 0.0;
		
};

 double CustomInitialCondition_rho_v_x::value(double x, double y) const 
	{     

		
		double a = calculate_a(x,y, this->kappa); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));
double q_2 = 2./(this->kappa-1.) * (1-a*a);
double q = Hermes::sqrt(q_2);

return q*rho;//Hermes::sqrt(2.);


			
};

 Ord CustomInitialCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_x::clone() const {     
			return new CustomInitialCondition_rho_v_x(this->mesh, kappa);
 };
//-------rho_v_y
 void CustomInitialCondition_rho_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
		dx =0.;
		dy = 0.0;
		
};

 double CustomInitialCondition_rho_v_y::value(double x, double y) const 
	{     

		
		double a = calculate_a(x,y, this->kappa); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));
double q_2 = 2./(this->kappa-1.) * (1-a*a);
double q = Hermes::sqrt(q_2);

//return -q*rho/Hermes::sqrt(2.);
return 0;

			
};

 Ord CustomInitialCondition_rho_v_y::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* CustomInitialCondition_rho_v_y::clone() const {     
			return new CustomInitialCondition_rho_v_y(this->mesh, kappa);
 };



//----------energy----------


 void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
	
		dx =0.;
		dy = 0.0;	
	
};

 double CustomInitialCondition_e::value(double x, double y) const 
{  
		
		double a = calculate_a(x,y, this->kappa); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));
double q_2 = 2./(this->kappa-1.) * (1-a*a);
double q = Hermes::sqrt(q_2);

double pressure = 1./this->kappa*std::pow(a, (2.*this->kappa/(this->kappa-1.))); 



double rho_v_x = rho*q/Hermes::sqrt(2.);

	return QuantityCalculator::calc_energy(rho, rho_v_x ,-rho_v_x, pressure, this->kappa);

};

 Ord CustomInitialCondition_e::ord(double x, double y)   const {
      return Ord(4);
};
 MeshFunction<double>* CustomInitialCondition_e::clone() const
    {
return new CustomInitialCondition_e(this->mesh,this->kappa);

    }

//------bdry_rho
 void BoundaryCondition_rho::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho::value(double x, double y) const 
	{    		


	double q = 0.5;
	double a = std::sqrt(1.-(this->kappa-1.)/2.*q*q); 
	double	rho = std::pow(a, (2./(this->kappa-1.)));

return rho;	


};

 Ord BoundaryCondition_rho::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho::clone() const {     
			return new BoundaryCondition_rho(this->mesh, this->kappa);
 };

// v_x
 void BoundaryCondition_rho_v_x::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho_v_x::value(double x, double y) const 
	{    			
		
	double q = 0.5;
	double a = std::sqrt(1.-(this->kappa-1.)/2.*q*q); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));

//return q*rho/Hermes::sqrt(2.);

return 0;
};

 Ord BoundaryCondition_rho_v_x::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_x::clone() const {     
			return new BoundaryCondition_rho_v_x(this->mesh, this->kappa);
 };
// v_y
 void BoundaryCondition_rho_v_y::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho_v_y::value(double x, double y) const 
	{    			
		
	double q = 0.5;
	double a = std::sqrt(1.-(this->kappa-1.)/2.*q*q); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));

//return -q*rho/Hermes::sqrt(2.);
return 0;
};

 Ord BoundaryCondition_rho_v_y::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_v_y::clone() const {     
			return new BoundaryCondition_rho_v_y(this->mesh, this->kappa);
 };


//bry_e
 void BoundaryCondition_rho_e::derivatives(double x, double y, double& dx, double& dy) const 
{ 
	
		dx =0.;
		dy = 0.0;	

};

 double BoundaryCondition_rho_e::value(double x, double y) const 
	{    			

	double q = 0.5;
	double a = std::sqrt(1.-(this->kappa-1.)/2.*q*q); 
		double		rho = std::pow(a, (2./(this->kappa-1.)));

double pressure = 1./this->kappa*std::pow(a, (2.*this->kappa/(this->kappa-1.))); 



double rho_v_x = rho*q/Hermes::sqrt(2.);

	return QuantityCalculator::calc_energy(rho, rho_v_x ,-rho_v_x, pressure, this->kappa);



};

 Ord BoundaryCondition_rho_e::ord(double x, double y)   const {
      return Ord(4);
};

 MeshFunction<double>* BoundaryCondition_rho_e::clone() const {     
			return new BoundaryCondition_rho_e(this->mesh, kappa);
 };

