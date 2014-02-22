
bool equal(double x, double y)
{
		if(std::fabs(x-y)< 1e-10) return true;
		else return false;
}

int boundary_test(double x, double y, double kappa)
{
	if((equal(y,0.))) return 4;
	double a, rho,q, p, J, k ;
	double vn_x[4], vn_y[4];
	q= 0.5;
	a = std::sqrt(1.-( kappa-1.)/2.*q*q);
	rho = std::pow(a, (2./( kappa-1.)));
	p = 1./ kappa*std::pow(a, (2.* kappa/( kappa-1.)));
	J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
	k = 1.5;
	vn_x[3] = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;				
	vn_y[3]= 1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));				
	k = 0.7;
	vn_x[2] = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;
	vn_y[2]= 1./(k*q*rho)*std::sqrt(1.-std::pow(q/k,2.));

	if(y<vn_y[3]) return 0.;
	if(x>vn_x[2]) return 0.;

	return 3;

}

double calculate_x(double k, double y, double kappa)
{	
	double q = 0.5;	
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
double rho, q_2, g_a,drho,dq_2,dg_a, a_old;
	for(int i = 0; i<100;i++)
	{
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
			
		 g_a = y*y-1./(k*k*rho*rho*q_2) +1./(std::pow(k,4.)*rho*rho) ;
		if(std::fabs(g_a)<1e-10) break;
		
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 dg_a= 	dq_2/(k*k*rho*rho*q_2*q_2)+2.*drho/(k*k*q_2*std::pow(rho,3.))-2.*drho/(std::pow(k,4.)*std::pow(rho,3.));

		 a_old = a;
		a = a_old - g_a/dg_a;

	}

		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
	double J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));	
double x = 1./(2.*rho)*(2./(k*k)-1./q_2)-J/2.;
return x; 
}


double calculate_y(double k, double x, double kappa)
{	
	double q = 0.5;	
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
	double rho, q_2, g_a,drho,dq_2,dg_a, a_old,J,dJ;
	for(int i = 0; i<100;i++)
	{
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
		 J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));		
		 g_a = x - 1./(2.*rho)*(2./(k*k)-1./q_2)+J/2.;
		if(std::fabs(g_a)<1e-10) break;
		 dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+std::log(std::exp(1))/(1-a*a));
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 dg_a= 	drho/(k*k*rho*rho)+0.5*(-drho/(rho*rho*q_2)-dq_2/(rho*q_2*q_2))+0.5*dJ;

		 a_old = a;
		a = a_old - g_a/dg_a;

	}

	rho = std::pow(a, (2./(kappa-1.)));
	q_2 = 2./(kappa-1.) * (1-a*a);
	q = std::sqrt(q_2);
	double y = 1./(k*rho*q) *std::sqrt(1-q_2/(k*k));
	return y; 
}

/*
void test_best(double x, double y, double k, double kappa, double& nx, double& ny)
{


	double q = (0.5+k)/2.;	
	double a =std::sqrt(1.-(kappa-1.)*0.5*q*q);
double k_2 = k*k;
	double dist = 0.1;
double x_a, y_a, dx_da, dy_da, dxx_da2, dyy_da2;
double rho, rho_2, q_2, drho,dq_2,J,dJ, dJJ, ddrho,ddq_2, dq;
double res_1, res_2;
double a_old, dist_old, q_old, res_old;


double a_11, a_12, a_21, a_22;
double b_1, b_2; double x_1,x_2;


for(int i = 0; i<100; i++)
{

		 q_2 = 2./(kappa-1.) * (1-a*a);
		 q = std::sqrt(q_2);
if(q>=k){ q = k;q_2 = q*q;  a =std::sqrt(1.-(kappa-1.)*0.5*q_2);};

	rho = std::pow(a, (2./(kappa-1.)));
	rho_2 = rho*rho; 
		 J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
	x_a = 0.5/rho*(2./k_2 - 1./q_2)-J/2.;
	if(q<k) y_a = 1./(k*rho*q)*std::sqrt(1.-q_2/k_2);
	else y_a = 0;

		 dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+1./(1-a*a));
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		dq = -2.*a/std::sqrt(2.*(kappa-1.)*(1-a*a));

	dx_da = -(drho/(k_2*rho_2)+0.5*(-drho/(rho_2*q_2)-dq_2/(rho*q_2*q_2))+0.5*dJ);
	if(q<k) dy_da = 0.5/std::sqrt(k_2-q_2)*(-2.*k_2*drho*q_2-k_2*rho*dq_2+2.*drho*q_2*q_2)/(k_2*rho_2*q*q_2);
	else dy_da = 0;


res_1 = x_a + dist* dy_da - x;
res_2 = y_a - dist* dx_da - y;
if( (res_1*res_2+res_2*res_2)<1e-8) 
{//printf("small enough:%e",(res_1*res_2+res_2*res_2) );
	break;}
//if(q>=k) { printf("q>= k :%e, res = %f, res_old = %f", q,(res_1*res_2+res_2*res_2), res_old); break;  }
//if(q<0.5) { printf("q<k :%e, res = %f, res_old = %f", q,(res_1*res_2+res_2*res_2), res_old);  break; }
if((q>=k)||(q<0.5)) break;


dJJ = 2.*(1./std::pow(a,3) +2./std::pow(a,5)+3./std::pow(a,7)-a/Hermes::sqr(1-a*a));
ddrho = 2.*(3.-kappa)/Hermes::sqr(kappa-1.)*std::pow(a, (4.-2.*kappa)/(kappa-1));
ddq_2 = -4./(kappa-1.);

dxx_da2=-(1./(k_2)*(ddrho/(rho_2) - 2.* drho*drho/(std::pow(rho,3.))) 
			+ 0.5*(-ddrho/(rho_2*q_2)+2.* drho*drho/(std::pow(rho,3.)*q_2)+ drho*dq_2/(rho_2*q_2*q_2))   
			+0.5*(-ddq_2/(rho*q_2*q_2)+ drho*dq_2/(rho_2*q_2*q_2)+ 2.* dq_2*dq_2/(rho*std::pow(q_2,3.))) +0.5*dJJ);


dyy_da2 = - dq_2*0.25/std::pow(k_2-q_2, 1.5) * (2.*k_2*drho*q_2+k_2*rho*dq_2-2.*drho*q_2*q_2)/(k_2*rho_2*q*q_2) 
			- 0.5/std::sqrt(k_2-q_2)*(2.*k_2*(ddrho*q_2+drho*dq_2)+k_2*(drho*dq_2+rho*ddq_2)-2.*(ddrho*q_2*q_2+drho*2.*q_2*dq_2))/(k_2*rho_2*q*q_2)
			- 0.5/std::sqrt(k_2-q_2)*(2.*k_2*drho*q_2+k_2*rho*dq_2-2.*drho*q_2*q_2)/k_2*(-2.*drho/std::pow(rho*q,3.)- 3*dq/(rho_2*q_2*q_2));

a_11= dx_da + dist*dyy_da2; a_12 = dy_da;
a_21 = dy_da - dist*dxx_da2; a_22 = -dx_da;
b_1 = -res_1; b_2= -res_2; 
//printf("jac = %f,%f,\n %f, %f \n", a_11,a_12,a_21,a_22 );
//printf("q= %f, d = %f \n", q, dist );
//solve via Gauss
double factor = a_21/a_11;
a_22-= factor*a_12;
b_2 -= factor*b_1;

x_2 = b_2/a_22;
x_1 = (b_1-a_12*x_2)/a_11;

a_old = a; dist_old = dist; q_old = q; res_old = (res_1*res_2+res_2*res_2);
dist+=x_2; a += x_1;




}

printf("q = %f, a= %f, d=  %f",q, a,dist);

double norm = std::sqrt(dy_da*dy_da+dx_da*dx_da);
if(q = k){
nx = 1.; ny = 0;
}else{
nx = dy_da/norm;
ny = -dx_da/norm;
}
if(x>0.1)
	if(nx<0){nx*=-1, ny*=-1;}

}

*/


double calculate_x_a(double k, double  kappa, double a)
{double k_2 = k*k;
	double	rho = std::pow(a, (2./(kappa-1.)));
double q_2 = 2./(kappa-1.) * (1-a*a);

	double 	J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
return (0.5/rho*(2./k_2 - 1./q_2)-J/2.);

}

double calculate_y_a(double k, double  kappa, double a)
{
double k_2 = k*k;
	double	rho = std::pow(a, (2./(kappa-1.)));
double q_2 = 2./(kappa-1.) * (1-a*a);
double	q = std::sqrt(q_2);
if(equal(k_2,q_2)) return 0;
return (1./(k*rho*q)*std::sqrt(1.-q_2/k_2));

}

double calculate_dx_da(double k, double  kappa, double a)
{
double rho = std::pow(a, (2./(kappa-1.)));
double k_2 = k*k;
double q_2 = 2./(kappa-1.) * (1-a*a);
double	q = std::sqrt(q_2);
if(equal(q,k)) return 1.;
double rho_2 = rho*rho; 
		double	dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+1./(1-a*a));
		double	drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		double	dq_2 = -4./(kappa-1.) * a;

return -(drho/(k_2*rho_2)+0.5*(-drho/(rho_2*q_2)-dq_2/(rho*q_2*q_2))+0.5*dJ);
}

double calculate_dy_da(double k, double  kappa, double a)
{

double rho = std::pow(a, (2./(kappa-1.)));
double k_2 = k*k;
double q_2 = 2./(kappa-1.) * (1-a*a);
double	q = std::sqrt(q_2);
if(equal(q,k)) return 0.;
double rho_2 = rho*rho; 
				double	drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		double	dq_2 = -4./(kappa-1.) * a;

return   0.5/std::sqrt(k_2-q_2)*(-2.*k_2*drho*q_2-k_2*rho*dq_2+2.*drho*q_2*q_2)/(k_2*rho_2*q*q_2);
}

double dx_da(double x, double k, double kappa)
{

//erst einmal a bestimmen
	double q = 0.5;	
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
double rho, q_2, g_a,drho,dq_2,dg_a, a_old,J,dJ;
	for(int i = 0; i<100;i++)
	{
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
		 J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));		
		 g_a = x - 1./(2.*rho)*(2./(k*k)-1./q_2)+J/2.;
		if(std::fabs(g_a)<1e-10) break;
		 dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+std::log(std::exp(1))/(1-a*a));
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 dg_a= 	drho/(k*k*rho*rho)+0.5*(-drho/(rho*rho*q_2)-dq_2/(rho*q_2*q_2))+0.5*dJ;

		 a_old = a;
		a = a_old - g_a/dg_a;
	}
// a bestimmt nun dx_da berechnen
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
		 dJ = -(1./(a*a) + 1./(std::pow(a,4)) +1./(std::pow(a,6))+std::log(std::exp(1))/(1-a*a));
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 dg_a= 	-(drho/(k*k*rho*rho)+0.5*(-drho/(rho*rho*q_2)-dq_2/(rho*q_2*q_2))+0.5*dJ);

return dg_a;
}

double dy_da(double y, double k, double kappa)
{
//erst einmal a bestimmen
	double q = 0.5;	
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
double rho, q_2, g_a,drho,dq_2,dg_a, a_old;
	for(int i = 0; i<100;i++)
	{
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
			
		 g_a = y*y-1./(k*k*rho*rho*q_2) +1./(std::pow(k,4.)*rho*rho) ;
		if(std::fabs(g_a)<1e-10) break;
		
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		 dg_a= 	dq_2/(k*k*rho*rho*q_2*q_2)+2.*drho/(k*k*q_2*std::pow(rho,3.))-2.*drho/(std::pow(k,4.)*std::pow(rho,3.));

		 a_old = a;
		a = a_old - g_a/dg_a;

	}
//dy/da bestimmen
		rho = std::pow(a, (2./(kappa-1.)));
		 q_2 = 2./(kappa-1.) * (1-a*a);
		q = Hermes::sqrt(q_2);
		 drho = 2./(kappa-1.)*std::pow(a, (2./(kappa-1.)-1.));
		 dq_2 = -4./(kappa-1.) * a;
		dg_a = 0.5*rho*q/(Hermes::sqrt(k*k-q_2))*(-2.*drho/(std::pow(rho,3)*q_2)-dq_2/(rho*rho*q_2*q_2)+2.*drho/(k*k*std::pow(rho,3)));

return dg_a; 
}



double dx_dk(double x, double q, double kappa)
{
//erst einmal k bestimmen:
double k = 1.0;
double g_k,dg_k, k_old;
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
	double	rho = std::pow(a, (2./(kappa-1.)));
	double	 q_2 = q*q;	
double J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));	

	for(int i = 0; i<100;i++)
	{		
		 g_k = x - (1./(rho*k*k)-0.5/(rho*q_2)-0.5*J);
		if(std::fabs(g_k)<1e-10) break;
		 dg_k= 	2./(rho*std::pow(k,3.));
		 k_old = k;
		k = k_old - g_k/dg_k;
	}
//dx/dk bestimmen
 dg_k =- 2./(rho*std::pow(k,3.));
return dg_k;

}

double dy_dk(double y, double q, double kappa)
{
//erst einmal k bestimmen:
double k = 1.0;
double g_k,dg_k, k_old;
	double a =std::sqrt(1.-(kappa-1.)/2.*q*q);
	double	rho = std::pow(a, (2./(kappa-1.)));
	double	 q_2 = q*q;	
double J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));	

	for(int i = 0; i<100;i++)
	{		
		 g_k = y - std::sqrt(k*k-q_2)/(k*k*rho*q);
		if(std::fabs(g_k)<1e-10) break;
		 dg_k= 	-(2.*q_2-k*k)/(std::sqrt(k*k-q_2)*std::pow(k,3.)*rho*q);
		 k_old = k;
		k = k_old - g_k/dg_k;
	}
//dx/dk bestimmen
 dg_k = (2.*q_2-k*k)/(std::sqrt(k*k-q_2)*std::pow(k,3.)*rho*q);
return dg_k;

}


