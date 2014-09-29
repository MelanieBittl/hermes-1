#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,int num_of_equations): WeakForm<double>(num_of_equations), num_of_equations(num_of_equations),prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy)
	{
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
		for(int k =0; k<4;k++)
		{
			add_vector_form(new EulerEquationsWeakForm_Mass::MassLinearform(k));
			add_matrix_form(new DefaultMatrixFormVol<double>(k,k)); 
		} 

	}


	    WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
    EulerEquationsWeakForm_Mass* wf;
    wf = new EulerEquationsWeakForm_Mass(this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy,4);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->ext[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        wf->ext.push_back(new Solution<double>());
        wf->ext.back()->copy(this->ext[i]);
      }
      else
        wf->ext.push_back(this->ext[i]->clone());
    }
    return wf;
    }
    
    

    double EulerEquationsWeakForm_Mass::MassLinearform::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		{
			double result = 0.;
          for (int k = 0; k < n; k++) {
            result += wt[k] * ext[this->i]->val[k]* v->val[k];
          }

		}


    Ord EulerEquationsWeakForm_Mass::MassLinearform::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
		{
		return Ord(10);
		}

    VectorFormVol<double>* EulerEquationsWeakForm_Mass::MassLinearform::clone() const
	{
	return new EulerEquationsWeakForm_Mass::MassLinearform(*this);
};
    


//_---------------Matrix K -------------------------------------------------------------------------------------------

	EulerK::EulerK(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy)
	{

	for(int k =0; k<4;k++)
	{	for(int i = 0; i<4;i++)			
			add_matrix_form(new EulerK::EulerEquationsBilinearForm(i,k,kappa));	
add_vector_form(new EulerK::EulerEquationsLinearForm(k,kappa));
}
    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

	};


	EulerK ::~EulerK ()
	{
		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerK::clone() const
    {
    EulerK* wf;
    wf = new EulerK(this->kappa, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->mirror_condition,4);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->ext[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        wf->ext.push_back(new Solution<double>());
        wf->ext.back()->copy(this->ext[i]);
      }
      else
        wf->ext.push_back(this->ext[i]->clone());
    }
    return wf;
    }


//---------------------K----------------

    double EulerK::EulerEquationsBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		 {
      double result = 0.;
      for (int i = 0;i < n;i++) 
      {

        result += wt[i] * u->val[i] *
        ( (static_cast<EulerK*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
          * v->dx[i]+        
         (static_cast<EulerK*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,entry_j) 
          * v->dy[i]);
		
      }
      return result;
    }
    


    Ord EulerK::EulerEquationsBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormVol<double>* EulerK::EulerEquationsBilinearForm::clone() const
    {
    return new EulerK::EulerEquationsBilinearForm(this->entry_i, this->entry_j, this->kappa);
    }



//------------K- Linearform
    double EulerK::EulerEquationsLinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
      double result = 0.;
      for (int i = 0;i < n;i++) 
      {

		for(int k =0; k<4;k++)
		{
		  result += wt[i] * ext[k]->val[i]*( v->dx[i]*
		   (static_cast<EulerK*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,k) 
		  + v->dy[i]*
		   (static_cast<EulerK*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,k)); 


		}
		
      }
      return result;
}


    Ord EulerK::EulerEquationsLinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(10);
    }

    VectorFormVol<double>* EulerK::EulerEquationsLinearForm::clone() const
    {
    return new EulerK::EulerEquationsLinearForm(this->entry_i, this->kappa);
    }


//---------------------------Boundary only------------------------
	EulerS::EulerS(double kappa,MeshSharedPtr mesh,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), rho_ext(rho_ext), v1_ext (v1_ext), v2_ext(v2_ext), energy_ext (energy_ext) , mesh(mesh)
	{

		for(int k =0; k<4;k++)
		{
			for(int i = 0; i<4;i++)
				add_matrix_form_surf(new EulerS::EulerBoundaryBilinearForm(kappa,i,k));
			
			add_vector_form_surf(new EulerS::EulerBoundaryLinearform(kappa, k));
		}

    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, rho_ext, v1_ext, v2_ext, energy_ext));

	};


	EulerS ::~EulerS ()
	{

		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerS::clone() const
    {
    EulerS* wf;
    wf = new EulerS(this->kappa,this->mesh,this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->mirror_condition,4);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->ext[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        wf->ext.push_back(new Solution<double>());
        wf->ext.back()->copy(this->ext[i]);
      }
      else
        wf->ext.push_back(this->ext[i]->clone());
    }
    return wf;
    }



bool equal_value(double x, double y)
{
		if(std::fabs(x-y)< 1e-10) return true;
		else return false;
}

int boundary(double x, double y, double kappa)
{
	if((equal_value(y,0.))) return 4;
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
/*	k=1.5;
	q = k;
	a = std::sqrt(1.-( kappa-1.)/2.*q*q);
	rho = std::pow(a, (2./( kappa-1.)));
	p = 1./ kappa*std::pow(a, (2.* kappa/( kappa-1.)));
	J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
	vn_x[0] = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;				
	vn_y[0]= 0;
	k = 0.7;
	q = k;
	a = std::sqrt(1.-( kappa-1.)/2.*q*q);
	rho = std::pow(a, (2./( kappa-1.)));
	p = 1./ kappa*std::pow(a, (2.* kappa/( kappa-1.)));
	J = 1./a + 1./(3.*std::pow(a,3)) +1./(5.*std::pow(a,5))-0.5*std::log((1.+a)/(1-a));
	vn_x[1] = 1./(2.*rho)*(2./(k*k)-1./(q*q))-J/2.;				
	vn_y[1]= 0;


	if((equal_value(y,0.))&&(x>vn_x[0])&&(x<vn_x[1])) return 4;
	if(y<=vn_y[3]) return 0.;
	if(x>=vn_x[2]) return 0.;*/
	if(y<vn_y[3]) return 0.;
	if(x>vn_x[2]) return 0.;

		return 3;

}

double normalized(double x, double y)
{

	return Hermes::sqrt(x*x+y*y);
}


bool out_of_range(double k, double a, double kappa)
{
double q = std::sqrt(2./(kappa-1.) * (1-a*a));
	if(q>k) return true;
	if(q<0.5) return true;
return false;
}


double calculate_a_start(double x, double y, double kappa)
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



bool test_best(double x, double y, double k, double kappa, double& nx, double& ny)
{
	double q = 0.5;	
	//double a =std::sqrt(1.-(kappa-1.)*0.5*q*q);
	double a_start = calculate_a_start(x,y,kappa);
	double a = a_start;
	double k_2 = k*k;
	double dist = 0.1;
	double x_a, y_a, dx_da, dy_da, dxx_da2, dyy_da2;
	double rho, rho_2, q_2, drho,dq_2,J,dJ, dJJ, ddrho,ddq_2, dq;
	double res_1, res_2;
	double a_old, dist_old, q_old, res_old;
	double a_11, a_12, a_21, a_22;
	double b_1, b_2; double x_1,x_2;

bool test_1 = true;

int i = 0;

	do
	{		q_2 = 2./(kappa-1.) * (1-a*a);
			q = std::sqrt(q_2);
			if(q>=k){ if(test_1){ q= 0.5; test_1 = false;}
						else q = k;
						//q= k-0.001;
			q_2 = q*q;  a =std::sqrt(1.-(kappa-1.)*0.5*q_2);};
			if(q<0.5){if(test_1){q = k-0.001;test_1 = false;}
						else q= 0.5;	
					
					q_2 = q*q;  a =std::sqrt(1.-(kappa-1.)*0.5*q_2);}

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
			else{ dy_da = 1.; dx_da = 0;}


			res_1 = x_a + dist* dy_da - x;
			res_2 = y_a - dist* dx_da - y;
			if(normalized(res_1,res_2)<1e-8) 
			{//printf("small enough:%e",(res_1*res_2+res_2*res_2) );
			break;}
			//if(q>=k) { printf("q>= k :%e, res = %f, res_old = %f", q,(res_1*res_2+res_2*res_2), res_old); break;  }
			//if(q<0.5) { printf("q<k :%e, res = %f, res_old = %f", q,(res_1*res_2+res_2*res_2), res_old);  break; }
			//if((q>=k)||(q==0.5)){/*printf("exceed");*/break;}

			dJJ = 2.*(1./std::pow(a,3) +2./std::pow(a,5)+3./std::pow(a,7)-a/Hermes::sqr(1.-a*a));
			ddrho = 2.*(3.-kappa)/Hermes::sqr(kappa-1.)*std::pow(a, (4.-2.*kappa)/(kappa-1.));
			ddq_2 = -4./(kappa-1.);

			dxx_da2=-(1./(k_2)*(ddrho/(rho_2) - 2.* drho*drho/(std::pow(rho,3.))) 
				+ 0.5*(-ddrho/(rho_2*q_2)+2.* drho*drho/(std::pow(rho,3.)*q_2)+ drho*dq_2/(rho_2*q_2*q_2))   
				+0.5*(-ddq_2/(rho*q_2*q_2)+ drho*dq_2/(rho_2*q_2*q_2)+ 2.* dq_2*dq_2/(rho*std::pow(q_2,3.))) +0.5*dJJ);


			if(q<k) dyy_da2 = - dq_2*0.25/std::pow(k_2-q_2, 1.5) * (2.*k_2*drho*q_2+k_2*rho*dq_2-2.*drho*q_2*q_2)/(k_2*rho_2*q*q_2) 
				- 0.5/std::sqrt(k_2-q_2)*(2.*k_2*(ddrho*q_2+drho*dq_2)+k_2*(drho*dq_2+rho*ddq_2)-2.*(ddrho*q_2*q_2+drho*2.*q_2*dq_2))/(k_2*rho_2*q*q_2)
				- 0.5/std::sqrt(k_2-q_2)*(2.*k_2*drho*q_2+k_2*rho*dq_2-2.*drho*q_2*q_2)/k_2*(-2.*drho/std::pow(rho*q,3.)- 3*dq/(rho_2*q_2*q_2));
			else dyy_da2 = 0;

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

			a_old = a; dist_old = dist; q_old = q; res_old = (res_1*res_1+res_2*res_2);

		double damping=1.;
		/*for(int l = 0; l<100;l++)
		{
			if(!out_of_range(k, a+x_1*damping,kappa))break;
			else damping-=0.01;
		}*/

		dist+=damping*x_2; a += damping*x_1;

		i++;
	}while(i<100);


	//printf("q = %f, a= %f, d=  %f",q, a,dist);

	double norm = normalized(dy_da,dx_da);
	if(q == k){
	nx = 1.; ny = 0.;
//nx = 0.; ny = -1.;
	}else{
	nx = dy_da/norm;
	ny = -dx_da/norm;
	}
	if(k == 0.7)
		if(nx<0){nx*=-1, ny*=-1;}
	if(k == 1.5)
		if(nx>0){nx*=-1, ny*=-1;}

if(i==100) return false;
else return true;
}




//---Boudary Bilinearforms---------
   double EulerS::EulerBoundaryBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
	int bdry =0;
	double* ghost_state = new double[4];for(int l=0;l<4;l++) ghost_state[l]=0.;
	double * A_n = new double[4];
	double* dudu_j= new double[4];
	bool solid = false;
	double constant = 1.;
  double result = 0.;
double nx_ghost, ny_ghost, tx_ghost, ty_ghost, t_1,t_2, norm;
double nx, ny, tx, ty;
  for (int i = 0;i < n;i++) 
  {	
	bdry = boundary(e->x[i], e->y[i], this->kappa);
	if((bdry==4)&&(QuantityCalculator::calc_mach(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i],kappa)>1) ) 
//if((bdry==4)&&(QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i],kappa)>1) ) 
			bdry = 1;

	if(bdry== 0){solid = true; }
	else solid = false;

//if(solid==false) continue;


nx_ghost = e->nx[i];
ny_ghost = e->ny[i];
tx_ghost = e->tx[i];
ty_ghost = e->ty[i];


if( bdry==0)
{ double k;
	if(e->x[i]>0.1)
		k = 0.7;
	else 
		k = 1.5;
	bool converged = test_best(e->x[i], e->y[i], k, this->kappa, nx_ghost, ny_ghost);
	if(!converged) 
	{
	nx_ghost = e->nx[i];
	ny_ghost = e->ny[i];

	}
	tx_ghost = -ny_ghost; ty_ghost = nx_ghost;


}
nx = e->nx[i];
ny = e->ny[i];
tx = e->tx[i];
ty = e->ty[i];

nx= nx_ghost;
ny=	ny_ghost;
tx=	tx_ghost;
ty=	ty_ghost; 



		 if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 


(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state(bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx_ghost,ny_ghost,tx_ghost,ty_ghost, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);

if(bdry==1) constant =1.;
else constant = 0.5;


			if(bdry!=1){
					Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx,ny , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, entry_i,A_n); //ite-Zeile A_n

				if(bdry!=2) 
				{
						(static_cast<EulerS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx,ny ,tx, ty, ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,entry_j,dudu_j);//j-te Spalte
				}
			}

				
        result += wt[i] * u->val[i] *v->val[i]* constant*
        ( (static_cast<EulerS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
          * nx +
         (static_cast<EulerS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,entry_j)  
          * ny);
					 	

			if (bdry==2)
	      		result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
			else if(bdry!=1){ 
			 	result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
				for(int k =0;k<4;k++)
				{result += wt[i]*v->val[i] *u->val[i]*0.5*(
					((static_cast<EulerS*>(wf))->euler_fluxes->A(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],0,entry_i,k) 
									  * nx+
					(static_cast<EulerS*>(wf))->euler_fluxes->A(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],1,entry_i,k) 
									  * ny -A_n[k])* dudu_j[k]);
				} 
			}
			
		 }else{
				if(entry_i==1)
				{
						if(entry_j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*nx*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);		
	
						}else if(entry_j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*nx;
						}else if(entry_j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*nx;
						}else if(entry_j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*nx;
						}
				}else if(entry_i==2)
				{
						if(entry_j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*ny*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);			
						}else if(entry_j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*ny;
						}else if(entry_j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*ny;
						}else if(entry_j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ny;
						}

				}
			}
    }
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;

      return (-result);
  }  




    Ord EulerS::EulerBoundaryBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerS::EulerBoundaryBilinearForm::clone() const
    {
					return new EulerS::EulerBoundaryBilinearForm(kappa,this->entry_i, this->entry_j);
    }





//------------------------------------------------------------------
//----------------------------Linearform Boundary---------------------------------
//----------------------------------------------


    double EulerS::EulerBoundaryLinearform::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		 {
		double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
		double lambda =0.;
		double* new_variables = new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
int bdry; bool solid = false;
  double result = 0.;
double nx_ghost, ny_ghost, tx_ghost, ty_ghost, t_1,t_2, norm, nx, ny;
  for (int i = 0;i < n;i++) 
{
	
bdry = boundary(e->x[i], e->y[i], this->kappa);
	if((bdry==4)&&(QuantityCalculator::calc_mach(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i],kappa)>1) ) 
//if((bdry==4)&&(QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i],kappa)>1) ) 
			bdry = 1;
	if(bdry== 0){ solid = true; }
	else solid = false;




nx_ghost = e->nx[i];
ny_ghost = e->ny[i];
tx_ghost = e->tx[i];
ty_ghost = e->ty[i];





if( bdry==0)
{ double k;
	if(e->x[i]>0.1)
		k = 0.7;
	else 
		k = 1.5;
	bool converged = test_best(e->x[i], e->y[i], k, this->kappa, nx_ghost, ny_ghost);
if(!converged) 
{
nx_ghost = e->nx[i];
ny_ghost = e->ny[i];

}

tx_ghost = -ny_ghost; ty_ghost = nx_ghost;

}
nx = e->nx[i];
ny = e->ny[i];

nx= nx_ghost ;
ny=ny_ghost;


			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];

    if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];			

(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx_ghost,ny_ghost,tx_ghost,ty_ghost, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], new_variables, solid);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];	

				for(int k =0; k<4;k++)
				{
							  result += wt[i] *nx*v->val[i]*0.5* (new_variables[k]*
							   (static_cast<EulerS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,entry_i,k) 
									+	ext[k]->val[i]
							  * (static_cast<EulerS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,k)) ;
							  result += wt[i]* ny*v->val[i]*0.5*( new_variables[k]*
							   (static_cast<EulerS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,entry_i,k) 
										+ext[k]->val[i]
							  * (static_cast<EulerS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,k)) ;

				}

		if(bdry!=1)		
			result -= wt[i]*v->val[i]* 0.5 * 
								Boundary_helpers::calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, nx, ny,  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, entry_i);


				}else{//solid wall ->no mirror
						if(entry_i==1)
							result += wt[i]*v->val[i]*nx*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
						else if(entry_i==2)
							result += wt[i]*v->val[i]*ny*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
				}
			
			
		}
		delete [] new_variables;
     return (-result);
   }      



    Ord EulerS::EulerBoundaryLinearform::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* EulerS::EulerBoundaryLinearform::clone() const
    {
					return new EulerS::EulerBoundaryLinearform(this->kappa, this->entry_i);
    }

