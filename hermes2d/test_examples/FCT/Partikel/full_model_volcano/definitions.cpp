#include "definitions.h"

double y_top = 20;
double x_left = -14;
double x_right = 17;
double y_inlet = 1.;
int bdry_vent =  3;
int bdry_left = 1;
int bdry_right =1;
int bdry_top =1;
int bdry_rest = 0;

 EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(int num_of_equations): WeakForm<double>(num_of_equations), num_of_equations(num_of_equations)
	{
		for(int k =0; k<num_of_equations;k++)			
			add_matrix_form(new DefaultMatrixFormVol<double>(k,k)); 	 

	}
	    
WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
		EulerEquationsWeakForm_Mass* wf;
		wf = new EulerEquationsWeakForm_Mass(8);

    return wf;
    }

    


//_---------------Matrix K -------------------------------------------------------------------------------------------

	EulerK::EulerK(double gamma, MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(gamma)),gamma(gamma), riemann_invariants(new RiemannInvariants(gamma)),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g),
prev_density_p(prev_density_p), prev_density_vel_x_p(prev_density_vel_x_p), prev_density_vel_y_p(prev_density_vel_y_p), prev_energy_p(prev_energy_p)
	{
//gas phase
	for(int k =0; k<4;k++)
	{	for(int i = 0; i<4;i++)			
			add_matrix_form(new EulerK::EulerEquationsBilinearForm(i,k,gamma));	
	 add_vector_form(new EulerK::EulerEquationsLinearForm(k,gamma));
	}
//particle phase
	for(int k =4; k<8;k++)
	{	for(int i = 4; i<8;i++)			
			add_matrix_form(new EulerK::EulerEquationsBilinearForm(i,k,gamma,true));	
	 add_vector_form(new EulerK::EulerEquationsLinearForm(k,gamma,true));
	}
    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p, prev_density_vel_x_p, prev_density_vel_y_p, prev_energy_p));

	};


	EulerK ::~EulerK ()
	{
		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerK::clone() const
    {
    EulerK* wf;
    wf = new EulerK(this->gamma, this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,
this->prev_density_p, this->prev_density_vel_x_p, this->prev_density_vel_y_p, this->prev_energy_p,8);

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
		if(particle)
		{
			  for (int i = 0;i < n;i++) 
			  {
				result += wt[i] * u->val[i] *
				((static_cast<EulerK*>(wf))->euler_fluxes->A_p(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i], ext[7]->val[i],0,(entry_i-4),(entry_j-4)) 
				  * v->dx[i]+        
				 (static_cast<EulerK*>(wf))->euler_fluxes->A_p(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i], ext[7]->val[i],1,(entry_i-4),(entry_j-4)) 
				  * v->dy[i]);		
			  }
		}else{
			  for (int i = 0;i < n;i++) 
			  {
				result += wt[i] * u->val[i] *
				((static_cast<EulerK*>(wf))->euler_fluxes->A_g(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
				  * v->dx[i]+        
				 (static_cast<EulerK*>(wf))->euler_fluxes->A_g(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,entry_j) 
				  * v->dy[i]);		
			  }
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
    return new EulerK::EulerEquationsBilinearForm(this->entry_i, this->entry_j, this->gamma, this->particle);
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
			if(particle)
			{	  result += wt[i] * ext[k+4]->val[i]*( v->dx[i]*
					   (static_cast<EulerK*>(wf))->euler_fluxes->A_p(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i], ext[7]->val[i],0,(entry_i-4),k) 
					  + v->dy[i]*
					   (static_cast<EulerK*>(wf))->euler_fluxes->A_p(ext[4]->val[i], ext[5]->val[i], ext[6]->val[i], ext[7]->val[i],1,(entry_i-4),k));	  

			}else{
					  result += wt[i] * ext[k]->val[i]*( v->dx[i]*
					   (static_cast<EulerK*>(wf))->euler_fluxes->A_g(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,k) 
					  + v->dy[i]*
					   (static_cast<EulerK*>(wf))->euler_fluxes->A_g(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,k));
			}      
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
    return new EulerK::EulerEquationsLinearForm(this->entry_i, this->gamma, this->particle);
    }

   


//_---------------SourceTerm Part-------------------------------------------------------------------------------------------
	EulerSource::EulerSource(double particle_density, double d, double c_vg,double c_vp, double c_pg, double Pr, double mu, MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p, int num_of_equations): WeakForm<double>(num_of_equations),particle_density(particle_density),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g),
prev_density_p(prev_density_p), prev_density_vel_x_p(prev_density_vel_x_p), prev_density_vel_y_p(prev_density_vel_y_p), prev_energy_p(prev_energy_p),
 d(d), c_vg(c_vg),c_vp(c_vp), c_pg(c_pg), Pr(Pr), mu(mu)
	{
//gas/particle phase
	for(int k =0; k<8;k++)
	{	for(int i = 0; i<8;i++)			
			add_matrix_form(new EulerSource::EulerSourceBilinearForm(i,k));	
	 add_vector_form(new EulerSource::EulerSourceLinearForm(k));
	}    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p, prev_density_vel_x_p, prev_density_vel_y_p, prev_energy_p));

	};


	EulerSource ::~EulerSource ()
	{
	};
	
	WeakForm<double>* EulerSource::clone() const
    {
    EulerSource* wf;
    wf = new EulerSource(this->particle_density,d, c_vg,c_vp,c_pg,Pr,mu, this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,
this->prev_density_p, this->prev_density_vel_x_p, this->prev_density_vel_y_p, this->prev_energy_p,8);

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


//------------------Jacobian source----------------

    double EulerSource::EulerSourceBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		 {
if((entry_i==0)||(entry_i==4)) return 0;
      double result = 0.;

		double density_particle= (static_cast<EulerSource*>(wf))->particle_density;
		double diameter = (static_cast<EulerSource*>(wf))->d;		
		double c_vg = (static_cast<EulerSource*>(wf))->c_vg;
		double c_vp = (static_cast<EulerSource*>(wf))->c_vp;
		double c_pg = (static_cast<EulerSource*>(wf))->c_pg;
		double Pr =(static_cast<EulerSource*>(wf))->Pr;		
		double mu = (static_cast<EulerSource*>(wf))->mu;
		double kap = c_pg*mu/Pr;
      
		for (int i = 0;i < n;i++){
			double alpha_p  =ext[4]->val[i]/density_particle ;
			double alpha_g = 1.-alpha_p;
			
				double u_g_1 = ext[0]->val[i];
						double u_g_2 = ext[1]->val[i];
						double u_g_3 = ext[2]->val[i];
						double u_g_4 = ext[3]->val[i];
						
						double u_p_1 = ext[4]->val[i];
						double u_p_2 = ext[5]->val[i];
						double u_p_3 = ext[6]->val[i];
						double u_p_4 = ext[7]->val[i];
		
						double rho_g = u_g_1/alpha_g;  
						double rho_v_x_g = u_g_2/alpha_g; 
						double rho_v_y_g = u_g_3/alpha_g; 
						double rho_e_g = u_g_4/alpha_g;
						double v_x_g = rho_v_x_g/rho_g;
						double v_y_g = rho_v_y_g/rho_g;

						double rho_p = density_particle;  
						double rho_v_x_p = u_p_2/alpha_p; 
						double rho_v_y_p = u_p_3/alpha_p; 
						double rho_e_p = u_p_4/alpha_p;
						double v_x_p = rho_v_x_p/rho_p;
						double v_y_p = rho_v_y_p/rho_p;

						double v1_diff = (v_x_g - v_x_p);
						double v2_diff = (v_y_g - v_y_p);
						double v_diff_abs = std::sqrt(v1_diff*v1_diff+ v2_diff*v2_diff);

						double Re = rho_g*diameter*v_diff_abs/mu;
						if(Re==1) printf("Reynolds gleich 0!!!!!");
						double C_D = 0.44;
						if(Re<1000)
								C_D=24./Re*(1.+0.15*std::pow(Re,0.687));
						double Nu = 2.+0.65*std::sqrt(Re)*std::pow(Pr,1./3.);

						double T_g = 1./c_vg*(rho_e_g/rho_g-0.5*(v_x_g*v_x_g+v_y_g*v_y_g));
						double T_p = 1./c_vp*(rho_e_p/rho_p-0.5*(v_x_p*v_x_p+v_y_p*v_y_p));
						
						double Q_drag = 0.75*v_diff_abs*C_D/(diameter*alpha_g*density_particle);
						double Q_tem = Nu*6.*kap/(diameter*diameter*density_particle);


			if(entry_i==1)
			{
				if(entry_j==0)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_p_2);
				else if(entry_j==1)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_1;
				else if(entry_j==4)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*u_g_2;
				else if(entry_j==5)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_g_1);

			}else if(entry_i==2)
			{	if(entry_j==0)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_p_3);
				else if(entry_j==2)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_1;
				else if(entry_j==4)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_3);
				else if(entry_j==6)
					result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_g_1);

			}else if(entry_i==3)
			{	if(entry_j==0)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(-(u_p_2*u_p_2+u_p_3*u_p_3))/u_p_1;
					result  -= wt[i] * u->val[i] *v->val[i] *Q_tem *(u_p_1/c_vg)*(-u_g_4/(u_g_1*u_g_1)+ (u_g_2*u_g_2+u_g_3*u_g_3)/(u_g_1*u_g_1*u_g_1));
				}else if(entry_j==1)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_2;
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*(u_p_1/c_vg)*(-u_g_2)/(u_g_1*u_g_1);
				}else if(entry_j==2)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_3;
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*(u_p_1/c_vg)*(-u_g_3)/(u_g_1*u_g_1);
				}else if(entry_j==3)
				{   result -= wt[i] * u->val[i] *v->val[i]*Q_tem/u_g_1*(u_p_1/c_vg);
				}else if(entry_j==4)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_1*(u_p_2*u_p_2+u_p_3*u_p_3)/(u_p_1*u_p_1));
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*(T_g-T_p+u_p_4/(c_vp*u_p_1)-(u_p_2*u_p_2+u_p_3*u_p_3)/(c_vp*u_p_1*u_p_1));
				}else if(entry_j==5)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_2*u_p_1-2.*u_g_1*u_p_2)/u_p_1;
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*u_p_2/(c_vp*u_p_1);
				}else if(entry_j==6)
				{	result  -= wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_3*u_p_1-2.*u_g_1*u_p_3)/u_p_1;
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*u_p_3/(c_vp*u_p_1);
				}else if(entry_j==7)
					result  -= wt[i] * u->val[i] *v->val[i]*Q_tem*(-1./c_vp);
			}else if(entry_i==5)
			{
				if(entry_j==0)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_p_2);
				else if(entry_j==1)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_1;
				else if(entry_j==4)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*u_g_2;
				else if(entry_j==5)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_g_1);
			}else if(entry_i==6)
			{	if(entry_j==0)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_p_3);
				else if(entry_j==2)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_1;
				else if(entry_j==4)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_3);
				else if(entry_j==6)
					result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(-u_g_1);
			}else if(entry_i==7)
			{	if(entry_j==0)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(-(u_p_2*u_p_2+u_p_3*u_p_3))/u_p_1;
					result  += wt[i] * u->val[i] *v->val[i] *Q_tem *(u_p_1/c_vg)*(-u_g_4/(u_g_1*u_g_1)+ (u_g_2*u_g_2+u_g_3*u_g_3)/(u_g_1*u_g_1*u_g_1));
				}else if(entry_j==1)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_2;
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*(u_p_1/c_vg)*(-u_g_2)/(u_g_1*u_g_1);
				}else if(entry_j==2)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*u_p_3;
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*(u_p_1/c_vg)*(-u_g_3)/(u_g_1*u_g_1);
				}else if(entry_j==3)
				{   result += wt[i] * u->val[i] *v->val[i]*Q_tem/u_g_1*(u_p_1/c_vg);
				}else if(entry_j==4)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_1*(u_p_2*u_p_2+u_p_3*u_p_3)/(u_p_1*u_p_1));
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*(T_g-T_p+u_p_4/(c_vp*u_p_1)-(u_p_2*u_p_2+u_p_3*u_p_3)/(c_vp*u_p_1*u_p_1));
				}else if(entry_j==5)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_2*u_p_1-2.*u_g_1*u_p_2)/u_p_1;
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*u_p_2/(c_vp*u_p_1);
				}else if(entry_j==6)
				{	result  += wt[i] * u->val[i] *v->val[i] *Q_drag*(u_g_3*u_p_1-2.*u_g_1*u_p_3)/u_p_1;
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*u_p_3/(c_vp*u_p_1);
				}else if(entry_j==7)
					result  += wt[i] * u->val[i] *v->val[i]*Q_tem*(-1./c_vp);
			}			

		}

      return (result);
    }
    


    Ord EulerSource::EulerSourceBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormVol<double>* EulerSource::EulerSourceBilinearForm::clone() const
    {
    return new EulerSource::EulerSourceBilinearForm(this->entry_i, this->entry_j);
    }



//------------source- Linearform
    double EulerSource::EulerSourceLinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{

	if((entry_i==0)||(entry_i==4)) return 0;
double density_particle = (static_cast<EulerSource*>(wf))->particle_density;
		double diameter = (static_cast<EulerSource*>(wf))->d;		
		double c_vg = (static_cast<EulerSource*>(wf))->c_vg;
		double c_vp = (static_cast<EulerSource*>(wf))->c_vp;
		double c_pg = (static_cast<EulerSource*>(wf))->c_pg;
		double Pr =(static_cast<EulerSource*>(wf))->Pr;		
		double mu = (static_cast<EulerSource*>(wf))->mu;
		double kap = c_pg*mu/Pr;
double F_D_1, F_D_2, Q_T;
      double result = 0.;
	for (int i = 0;i < n;i++){
		double alpha_p  =ext[4]->val[i]/density_particle ;
		double alpha_g = 1.-alpha_p;

						double u_g_1 = ext[0]->val[i];
						double u_g_2 = ext[1]->val[i];
						double u_g_3 = ext[2]->val[i];
						double u_g_4 = ext[3]->val[i];
						
						double u_p_1 = ext[4]->val[i];
						double u_p_2 = ext[5]->val[i];
						double u_p_3 = ext[6]->val[i];
						double u_p_4 = ext[7]->val[i];
		
						double rho_g = u_g_1/alpha_g;  
						double rho_v_x_g = u_g_2/alpha_g; 
						double rho_v_y_g = u_g_3/alpha_g; 
						double rho_e_g = u_g_4/alpha_g;
						double v_x_g = rho_v_x_g/rho_g;
						double v_y_g = rho_v_y_g/rho_g;

						double rho_p = density_particle;  
						double rho_v_x_p = u_p_2/alpha_p; 
						double rho_v_y_p = u_p_3/alpha_p; 
						double rho_e_p = u_p_4/alpha_p;
						double v_x_p = rho_v_x_p/rho_p;
						double v_y_p = rho_v_y_p/rho_p;

						double v1_diff = (v_x_g - v_x_p);
						double v2_diff = (v_y_g - v_y_p);
						double v_diff_abs = std::sqrt(v1_diff*v1_diff+ v2_diff*v2_diff);

						double Re = rho_g*diameter*v_diff_abs/mu;
						if(Re==1) printf("Reynolds gleich 0!!!!!");
						double C_D = 0.44;
						if(Re<1000)
								C_D=24./Re*(1.+0.15*std::pow(Re,0.687));
						double Nu = 2.+0.65*std::sqrt(Re)*std::pow(Pr,1./3.);

						double T_g = 1./c_vg*(rho_e_g/rho_g-0.5*(v_x_g*v_x_g+v_y_g*v_y_g));
						double T_p = 1./c_vp*(rho_e_p/rho_p-0.5*(v_x_p*v_x_p+v_y_p*v_y_p));
						
						double Q_drag = 0.75*v_diff_abs*C_D/(diameter*alpha_g*density_particle);
						double Q_tem = Nu*6.*kap/(diameter*diameter*density_particle);
						
						F_D_1 = Q_drag* (u_g_2*u_p_1-u_p_2*u_g_1);
						F_D_2 = Q_drag*(u_g_3*u_p_1-u_p_3*u_g_1); 
						Q_T = Q_tem*u_p_1*(T_g-T_p); 
						


//gas phase
			if(entry_i==1)  //vx			  
					 result -= wt[i]*v->val[i]*F_D_1;
			else if(entry_i==2)//vy
					result -= wt[i]*v->val[i]*F_D_2;
			else if(entry_i ==3)	//e			
				result -= wt[i]*v->val[i]*((F_D_1*v_x_p+F_D_2*v_y_p)+Q_T);	

//particle phase

			else if(entry_i==5)//vx			  
					  result += wt[i]*v->val[i]*F_D_1;
			else if(entry_i==6)//vy
					result += wt[i]*v->val[i]*F_D_2;
			else if(entry_i ==7)	//e			
				result += wt[i]*v->val[i]*((F_D_1*v_x_p+F_D_2*v_y_p)+Q_T);

		

				      
		}	
	      return result;
}


    Ord EulerSource::EulerSourceLinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(10);
    }

    VectorFormVol<double>* EulerSource::EulerSourceLinearForm::clone() const
    {
    return new EulerSource::EulerSourceLinearForm(this->entry_i);
    }





//---------------------------Boundary+Penalty ------------------------
	EulerBoundary::EulerBoundary(double gamma,double sigma,double particle_density,MeshFunctionSharedPtr<double>  rho_ext_g, MeshFunctionSharedPtr<double>  v1_ext_g, MeshFunctionSharedPtr<double>  v2_ext_g, MeshFunctionSharedPtr<double>  energy_ext_g,MeshFunctionSharedPtr<double>  rho_ext_p, MeshFunctionSharedPtr<double>  v1_ext_p, MeshFunctionSharedPtr<double>  v2_ext_p, MeshFunctionSharedPtr<double>  energy_ext_p, 
MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p,double eps,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(gamma)),gamma(gamma), riemann_invariants(new RiemannInvariants(gamma)),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g), 
	rho_ext_g(rho_ext_g), v1_ext_g(v1_ext_g), v2_ext_g(v2_ext_g), energy_ext_g(energy_ext_g),
   prev_density_p(prev_density_p), prev_density_vel_x_p(prev_density_vel_x_p), prev_density_vel_y_p(prev_density_vel_y_p), prev_energy_p(prev_energy_p), 
	rho_ext_p(rho_ext_p), v1_ext_p(v1_ext_p), v2_ext_p(v2_ext_p), energy_ext_p(energy_ext_p) ,sigma(sigma), eps(eps), 
	particle_density(particle_density)
	{
//gas phase
		for(int k =0; k<4;k++)
		{
			for(int i = 0; i<4;i++)
			{
				add_matrix_form_surf(new EulerBoundary::EulerBoundaryBilinearForm(gamma,i,k));
				if(((i==1) ||(i==2))&&((k==1) ||(k==2)))
					add_matrix_form_surf(new EulerBoundary::PenaltyBilinearForm(sigma,eps,i,k));	
			}				
			add_vector_form_surf(new EulerBoundary::EulerBoundaryLinearform(gamma, k));
			if((k==1) ||(k==2))
				add_vector_form_surf(new EulerBoundary::PenaltyLinearForm(sigma, k));	
		}

//particle phase
		for(int k =4; k<8;k++)
		{
			for(int i = 4; i<8;i++)
			{
				add_matrix_form_surf(new EulerBoundary::EulerBoundaryBilinearForm(gamma,i,k,true));
				if(((i==5) ||(i==6))&&((k==5) ||(k==6)))
					add_matrix_form_surf(new EulerBoundary::PenaltyBilinearForm(sigma,eps,i,k,true));	
			}

			add_vector_form_surf(new EulerBoundary::EulerBoundaryLinearform(gamma, k,true));
			if((k==5) ||(k==6))
				add_vector_form_surf(new EulerBoundary::PenaltyLinearForm(sigma, k,true));
		}

    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g, rho_ext_g, v1_ext_g, v2_ext_g, energy_ext_g, prev_density_p, prev_density_vel_x_p, prev_density_vel_y_p, prev_energy_p, rho_ext_p, v1_ext_p, v2_ext_p, energy_ext_p));

	};


	EulerBoundary ::~EulerBoundary ()
	{

		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerBoundary::clone() const
    {
		EulerBoundary* wf;
		wf = new EulerBoundary(this->gamma,this->sigma,this->particle_density,this->rho_ext_g, this->v1_ext_g, this->v2_ext_g, this->energy_ext_g,
			this->rho_ext_p, this->v1_ext_p, this->v2_ext_p, this->energy_ext_p,
	 this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,
	this->prev_density_p, this->prev_density_vel_x_p, this->prev_density_vel_y_p, this->prev_energy_p,8);

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


//-------Boudary Bilinearforms---------
   double EulerBoundary::EulerBoundaryBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
	int bdry =0;
	double* ghost_state = new double[4];for(int l=0;l<4;l++) ghost_state[l]=0.;
	double * A_n = new double[4];
	double* dudu_j= new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
	
	double constant = 1.;
  double result = 0.;
  for (int i = 0;i < n;i++) 
  {		
		bdry = bdry_rest;
		if(e->x[i]== x_left){
				bdry = bdry_left;
		}else if(e->y[i]==y_top){
			bdry = bdry_top;
		}else if (e->x[i]==x_right){
				bdry =bdry_right;			
		}else if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	
			bdry = bdry_vent;
	
		
					rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];
    
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];
		
		
		if(bdry!=0) 
			bdry = (static_cast<EulerBoundary*>(wf))->riemann_invariants->get_bdry_info_short( rho, rho_v_x, rho_v_y, rho_energy, e->nx[i],e->ny[i], e->tx[i], e->ty[i],
					rho_ext,rho_v_x_ext, rho_v_y_ext, rho_energy_ext,false);
			
		if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	
			bdry = bdry_vent;
		


	if(bdry==1) constant =1.;
	else constant = 0.5;

	//----------particle---------------------
	if(particle)
		{	if(bdry==0) continue;

			rho = ext[8]->val[i];  
			rho_v_x = ext[9]->val[i]; 
			rho_v_y = ext[10]->val[i]; 
			rho_energy = ext[11]->val[i];
    
				rho_ext = ext[12]->val[i];
				rho_v_x_ext = ext[13]->val[i];
				rho_v_y_ext = ext[14]->val[i];
				rho_energy_ext = ext[15]->val[i];	

			if((bdry ==3)||(bdry ==2)) //inlet
			{
				if(entry_i==entry_j)
				{
					(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state_p(bdry, rho, rho_v_x, rho_v_y, rho_energy, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, ghost_state);
					Boundary_helpers::calculate_A_p_n( rho, rho_v_x, rho_v_y, rho_energy, e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], gamma, (entry_i-4),A_n);				
					result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[(entry_j-4)];
				}
			}else constant = 1.;			

		result += wt[i] * u->val[i] *v->val[i]* constant*
		    ( (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho, rho_v_x, rho_v_y, rho_energy,0,(entry_i-4),(entry_j-4)) 
		      * e->nx[i] +
		     (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho, rho_v_x, rho_v_y, rho_energy,1,(entry_i-4),(entry_j-4))  
		      * e->ny[i]);

	}else{//-------------gas phase-------------------
			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];
    
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];

	 		if(bdry!=0){ 
				(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state(bdry,rho, rho_v_x, rho_v_y, rho_energy, e->nx[i],e->ny[i],e->tx[i],e->ty[i], rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, ghost_state);

				if(bdry!=1){
						Boundary_helpers::calculate_A_n(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], gamma, entry_i,A_n); //ite-Zeile A_n

					if(bdry!=2) 
					{
							(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_du_du(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,entry_j,dudu_j);//j-te Spalte
					}
				}

				
		    result += wt[i] * u->val[i] *v->val[i]* constant*
		    ( (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,entry_j) 
		      * e->nx[i] +
		     (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,entry_j)  
		      * e->ny[i]);
						 	

				if (bdry==2)
			  		result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
				else if(bdry!=1){ 
				 	result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
					for(int k =0;k<4;k++)
					{result += wt[i]*v->val[i] *u->val[i]*0.5*(
						((static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],0,entry_i,k) 
										  * e->nx[i]+
						(static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],1,entry_i,k) 
										  * e->ny[i] -A_n[k])* dudu_j[k]);
					} 
				}
			
			}else{ 
					if(entry_i==1)
					{
							if(entry_j==0){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*0.5*e->nx[i]*
							(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);	
	
							}else if(entry_j==1){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[1]->val[i]/ext[0]->val[i]*e->nx[i];
							}else if(entry_j==2){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[2]->val[i]/ext[0]->val[i]*e->nx[i];
							}else if(entry_j==3){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*e->nx[i];
							}
					}else if(entry_i==2)
					{
							if(entry_j==0){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*0.5*e->ny[i]*
							(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);			
							}else if(entry_j==1){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[1]->val[i]/ext[0]->val[i]*e->ny[i];
							}else if(entry_j==2){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[2]->val[i]/ext[0]->val[i]*e->ny[i];
							}else if(entry_j==3){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*e->ny[i];
							}

					}
			}
		}
    }
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;

      return (-result);
  }  




    Ord EulerBoundary::EulerBoundaryBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerBoundary::EulerBoundaryBilinearForm::clone() const
    {
					return new EulerBoundary::EulerBoundaryBilinearForm(gamma,this->entry_i, this->entry_j, this->particle);
    }





//------------------------------------------------------------------
//----------------------------Linearform Boundary---------------------------------
//----------------------------------------------


    double EulerBoundary::EulerBoundaryLinearform::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
 {
		double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
		double lambda =0.;
		double* ghost_state = new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
int bdry; 
  double result = 0.;
  for (int i = 0;i < n;i++) 
	{

		bdry = bdry_rest;
		if(e->x[i]== x_left){
				bdry = bdry_left;
		}else if(e->y[i]==y_top){
			bdry = bdry_top;
		}else if (e->x[i]==x_right){
				bdry =bdry_right;			
		}else if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	
			bdry = bdry_vent;
		
		
		
				if(bdry!=0) 
			bdry=(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_bdry_info_short( rho, rho_v_x, rho_v_y, rho_energy, 
										e->nx[i],e->ny[i], e->tx[i], e->ty[i],
					rho_ext,rho_v_x_ext, rho_v_y_ext, rho_energy_ext,false);
			
		if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	
			bdry = bdry_vent;

	//----------particle---------------------
	if(particle)
	{	
if(bdry==0) continue;	

			rho = ext[8]->val[i];  
			rho_v_x = ext[9]->val[i]; 
			rho_v_y = ext[10]->val[i]; 
			rho_energy = ext[11]->val[i];

    
				rho_ext = ext[12]->val[i];
				rho_v_x_ext = ext[13]->val[i];
				rho_v_y_ext = ext[14]->val[i];
				rho_energy_ext = ext[15]->val[i];	

(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state_p( bdry,rho,rho_v_x,rho_v_y,rho_energy,rho_ext,rho_v_x_ext,rho_v_y_ext,rho_energy_ext, ghost_state);

		rho_new=ghost_state[0]; 
		rho_v_x_new=ghost_state[1]; 
		rho_v_y_new=ghost_state[2]; 
		rho_energy_new=ghost_state[3];	
				for(int k =0; k<4;k++)
				{
				  result += wt[i] *e->nx[i]*v->val[i]*0.5* (ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,(entry_i-4),k) 
						+	ext[k+8]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho, rho_v_x, rho_v_y, rho_energy,0,(entry_i-4),k)) ;
				  result += wt[i]* e->ny[i]*v->val[i]*0.5*( ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,(entry_i-4),k) 
							+ext[k+8]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_p(rho, rho_v_x, rho_v_y, rho_energy,1,(entry_i-4),k)) ;

				}

				if((bdry ==3)||(bdry ==2)) //inlet	
					result -= wt[i]*v->val[i]* 0.5 * 
										Boundary_helpers::calculate_A_p_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, gamma, (entry_i-4));




	}else{//------------gas-------------
			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];

    
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];

		if(bdry!=0){ 
			(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state);

		rho_new=ghost_state[0]; 
		rho_v_x_new=ghost_state[1]; 
		rho_v_y_new=ghost_state[2]; 
		rho_energy_new=ghost_state[3];	
				for(int k =0; k<4;k++)
				{
				  result += wt[i] *e->nx[i]*v->val[i]*0.5* (ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,entry_i,k) 
						+	ext[k]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,k)) ;
				  result += wt[i]* e->ny[i]*v->val[i]*0.5*( ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,entry_i,k) 
							+ext[k]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,k)) ;

				}

				if(bdry!=1)		
					result -= wt[i]*v->val[i]* 0.5 * 
										Boundary_helpers::calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, gamma, entry_i);
		}else{//solid wall ->no mirror
				if(entry_i==1)
					result += wt[i]*v->val[i]*e->nx[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
				else if(entry_i==2)
					result += wt[i]*v->val[i]*e->ny[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
		}
	  }			
			
	 }
	delete [] ghost_state;
     return (-result);
  }      



    Ord EulerBoundary::EulerBoundaryLinearform::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* EulerBoundary::EulerBoundaryLinearform::clone() const
    {
					return new EulerBoundary::EulerBoundaryLinearform(this->gamma, this->entry_i, this->particle);
    }



//-------Penalty Bilinearforms---------
   double EulerBoundary::PenaltyBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.; double bdry =0;
if((entry_i==0) ||(entry_i==3)||(entry_i==4) ||(entry_i==7)) return 0.; 
if((entry_j==0) ||(entry_j==3)||(entry_j==4) ||(entry_j==7)) return 0.;
double nx, ny;

  for (int i = 0;i < n;i++) 
  {		
			nx = e->nx[i];
			ny = e->ny[i];

		
		if(e->x[i]== x_left){
				continue;
		}else if((e->y[i]==y_top)||(e->x[i]==x_right)){
				continue;		
		}else if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	continue;
		

	
		if(particle){
			

			double vn = (ext[9]->val[i]*nx + ext[10]->val[i]*ny)/ext[8]->val[i];
				double factor = ((2.*vn*vn+eps)/(std::sqrt(vn*vn+eps)));

				if(entry_i==5)
				{
						if(entry_j==5){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(nx)*factor;
						}else if(entry_j==6){
							result += wt[i]*v->val[i] *u->val[i]*ny*nx*factor;
						}
				}else if(entry_i==6)
				{
					 if(entry_j==5){
							result += wt[i]*v->val[i] *u->val[i]*ny*nx*factor;
						}else if(entry_j==6){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(ny)*factor;
						}
				}
		}else{
		double vn = (ext[1]->val[i]*nx+ ext[2]->val[i]*ny)/ext[0]->val[i];
		double factor = ((2.*vn*vn+eps)/(std::sqrt(vn*vn+eps)));

				if(entry_i==1)
				{
						if(entry_j==1){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(nx)*factor;
						}else if(entry_j==2){
							result += wt[i]*v->val[i] *u->val[i]*ny*nx*factor;
						}
				}else if(entry_i==2)
				{
					 if(entry_j==1){
							result += wt[i]*v->val[i] *u->val[i]*ny*nx*factor;
						}else if(entry_j==2){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(ny)*factor;
						}
				}
		}
			
    }


      return sigma*(-result);
  }  




    Ord EulerBoundary::PenaltyBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerBoundary::PenaltyBilinearForm::clone() const
    {
					return new EulerBoundary::PenaltyBilinearForm(sigma,eps,this->entry_i, this->entry_j, this->particle);
    }





//------------------------------------------------------------------
//----------------------------Linearform Penalty---------------------------------
//----------------------------------------------


    double EulerBoundary::PenaltyLinearForm::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
 {
if((entry_i==0) ||(entry_i==3)||(entry_i==4) ||(entry_i==7)) return 0.;
double bdry =0;
double material_density = (static_cast<EulerBoundary*>(wf))->particle_density;
  double result = 0.;
  for (int i = 0;i < n;i++) 
	{
			
		if(e->x[i]== x_left){
				continue;
		}else if((e->y[i]==y_top)||(e->x[i]==x_right)){
				continue;		
		}else if((e->y[i]==y_inlet)&&(e->x[i]>=1.2)&&(e->x[i]<=1.8))	continue;
 
		if(particle)
		{
			double vn = (ext[9]->val[i]*e->nx[i] + ext[10]->val[i]*e->ny[i]); 
			
				if(entry_i==5)
					result += wt[i]*v->val[i]*e->nx[i]*vn*std::fabs(vn);
				else if(entry_i==6)
					result += wt[i]*v->val[i]*e->ny[i]*vn*std::fabs(vn);

		}else{
			double vn = (ext[1]->val[i]*e->nx[i] + ext[2]->val[i]*e->ny[i]);	
				if(entry_i==1)
					result += wt[i]*v->val[i]*e->nx[i]*vn*std::fabs(vn);
				else if(entry_i==2)
					result += wt[i]*v->val[i]*e->ny[i]*vn*std::fabs(vn);
		}			
	 }

     return sigma*(-result);
  }      



    Ord EulerBoundary::PenaltyLinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* EulerBoundary::PenaltyLinearForm::clone() const
    {
					return new EulerBoundary::PenaltyLinearForm(this->sigma, this->entry_i, this->particle);
    }
