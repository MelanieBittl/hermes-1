#include "boundary_condition.h"

double calculate_lambda_max(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
					double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa){

	double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));

double lambda_1,lambda_2, lambda_3;
lambda_1 = fabs(n_x*v_x_mean + n_y*v_y_mean-c_mean);
lambda_2 = fabs(n_x*v_x_mean + n_y*v_y_mean);
lambda_3 = fabs(n_x*v_x_mean + n_y*v_y_mean+c_mean);

if(lambda_1>lambda_2)
		if(lambda_1>lambda_3) return lambda_1;
		else return lambda_3;
else 
		if(lambda_2>lambda_3) return lambda_2;
		else return lambda_3;
}


double calculate_A_n_U(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
					double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry){

	double u[4] = {rho_new-rho, rho_v_x_new-rho_v_x, rho_v_y_new-rho_v_y, rho_energy_new-rho_energy };

	if((u[0]==0.)&&(u[1]==0.)&&(u[2]==0.)&&(u[3]==0.)) return 0.;

	double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
	double v_n = v_x_mean*n_x+v_y_mean*n_y;
	double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
	double b_2 = (kappa-1)/(c_mean*c_mean);
	double b_1 = b_2*q;
	double lambda[4] ={fabs(v_n- c_mean),fabs(v_n),fabs(v_n + c_mean),fabs(v_n)};

//Eintraege siehe Diss Moeller Appendix C
	double R[4][4] = { 1, 1, 1, 0, 
										v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
										v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
										H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

	double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
										1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
										0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
										n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

	//lambda*L
//	for(int i =0;i<4;i++)
		//	for(int j=0;j<4;j++)
				//	L[i][j] *= lambda[i];			

	double result =0.;

	for(int i =0;i<4;i++)
			for(int j=0;j<4;j++)
					result +=R[entry][j]*L[j][i]*u[i]*lambda[j];

	return result;		


}


void calculate_A_n(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
					double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry_i, double* A_n){
	double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
	double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
	double v_n = v_x_mean*n_x+v_y_mean*n_y;
	double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
	double b_2 = (kappa-1)/(c_mean*c_mean);
	double b_1 = b_2*q;
	double lambda[4] ={fabs(v_n- c_mean),fabs(v_n),fabs(v_n + c_mean),fabs(v_n)};



//Eintraege siehe Diss Moeller Appendix C
	double R[4][4] = { 1, 1, 1, 0, 
										v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
										v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
										H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

	double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
										1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
										0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
										n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

for(int i =0;i<4;i++) A_n[i]=0.;

	//double result =0.;

	for(int i =0;i<4;i++)
			for(int j=0;j<4;j++){
					A_n[0] +=R[entry_i][j]*L[j][0]*lambda[j];
					A_n[1] +=R[entry_i][j]*L[j][1]*lambda[j];
					A_n[2] +=R[entry_i][j]*L[j][2]*lambda[j];
					A_n[3] +=R[entry_i][j]*L[j][3]*lambda[j];
				}

}

//-------Boundary Condition for Matrix S

	EulerBoundary::EulerBoundary(double kappa,Hermes::vector<MeshFunctionSharedPtr<double> >slns, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition), slns(slns), kappa(kappa)
	{

 		add_matrix_form_surf(new EulerBoundaryBilinearForm_rho(kappa,0,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
 		add_matrix_form_surf(new EulerBoundaryBilinearForm_rho(kappa,1,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
 		add_matrix_form_surf(new EulerBoundaryBilinearForm_rho(kappa,2,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
 		add_matrix_form_surf(new EulerBoundaryBilinearForm_rho(kappa,3,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));

    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_x( kappa,0,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_x( kappa,1,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_x( kappa,2,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_x( kappa,3,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));

    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_y( kappa,0,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_y( kappa,1,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_y( kappa,2,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_vel_y( kappa,3,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));

    add_matrix_form_surf(new EulerBoundaryBilinearForm_e( kappa,0,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_e( kappa,1,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_e( kappa,2,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_matrix_form_surf(new EulerBoundaryBilinearForm_e( kappa,3,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));

		add_vector_form_surf(new EulerBoundary_rho(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_vector_form_surf(new EulerBoundary_v_x( kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_vector_form_surf(new EulerBoundary_v_y( kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    add_vector_form_surf(new EulerBoundary_e(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition));
    
//this->slns = Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, rho_ext, v1_ext, v2_ext, energy_ext)

    this->set_ext(this->slns);

	};


	EulerBoundary ::~EulerBoundary ()
	{
		/*for(int i=0; i<this->vfsurf.size();i++)		
			delete get_vfsurf()[i];			
			for(int i=0; i<this->mfsurf.size();i++)
			delete get_mfsurf()[i];

		WeakForm<double>::delete_all();*/
		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	    WeakForm<double>* EulerBoundary::clone() const
    {
     // const_cast<EulerBoundary*>(this)->warned_nonOverride = false;
      //return new EulerBoundary(*this);

    EulerBoundary* wf;
    wf = new EulerBoundary(this->kappa, this->slns, this->mirror_condition);

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

    wf->set_current_time_step(this->get_current_time_step());

    return wf;


    }
	

//Bilinearforms
//rho

  template<typename Real, typename Scalar>
    Scalar EulerBoundaryBilinearForm_rho::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
{
	int bdry =0;
	double* ghost_state = new double[4];
	double * A_n = new double[4];
	double* dudu_j= new double[4];
		bool solid = true;
		double constant =0.5;
  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++) 
  {		

		 if((mirror_condition==true)||(solid==false)){ 
				bdry =riemann_invariants->get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


			if(bdry!=1){
						calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 0, A_n); //0te-Zeile A_n
				if(bdry!=2) 
						riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
			}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}

				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				(euler_fluxes->A_1_0_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_0_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				(euler_fluxes->A_1_0_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_0_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				(euler_fluxes->A_1_0_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_0_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				(euler_fluxes->A_1_0_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_0_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );

				}
			
		 }
    }
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;

      return (-result);
  }  


    double EulerBoundaryBilinearForm_rho::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerBoundaryBilinearForm_rho::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerBoundaryBilinearForm_rho::clone() const
    {
					return new EulerBoundaryBilinearForm_rho(kappa, this->j,this->euler_fluxes, this->riemann_invariants, this->mirror_condition);
    }

//vel_x
    template<typename Real, typename Scalar>
    Scalar EulerBoundaryBilinearForm_vel_x::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
{
	int bdry =0;
	double* ghost_state = new double[4];
	double * A_n = new double[4];
	double* dudu_j= new double[4];
		bool solid = true;
		double constant =0.5;
  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++) 
  {		
			if((mirror_condition==true)||(solid==false)){ 

	bdry =riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


	if(bdry!=1){
				calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 1, A_n); //1te-Zeile A_n
		if(bdry!=2) 
				riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
	}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] *constant
        * euler_fluxes->A_1_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				(euler_fluxes->A_1_1_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_1_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				(euler_fluxes->A_1_1_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_1_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				(euler_fluxes->A_1_1_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_1_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				(euler_fluxes->A_1_1_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_1_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );
				}
			}else{
						if(j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*e->nx[i]*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);		
	
						}else if(j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*e->nx[i];
						}else if(j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*e->nx[i];
						}else if(j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*e->nx[i];
						}
			}
    
		}

		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;
    return (-result);
 }

    double EulerBoundaryBilinearForm_vel_x::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerBoundaryBilinearForm_vel_x::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerBoundaryBilinearForm_vel_x::clone() const { return new EulerBoundaryBilinearForm_vel_x(kappa, this->j,this->euler_fluxes, this->riemann_invariants, this->mirror_condition); }

//vel_y
  template<typename Real, typename Scalar>
    Scalar EulerBoundaryBilinearForm_vel_y::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
		int bdry =0;
		double* ghost_state = new double[4];
		double* dudu_j = new double[4];
		double* A_n = new double[4];
		bool solid = true;
		double constant =0.5;

    Scalar result = Scalar(0);
    for (int i = 0;i < n;i++) 
    {				


			if((mirror_condition==true)||(solid==false)){ 
					bdry =riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);



					if(bdry!=1){
								calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 2, A_n); //2te-Zeile A_n
						if(bdry!=2) 
								riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
					}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				 if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				(euler_fluxes->A_1_2_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_2_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				(euler_fluxes->A_1_2_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_2_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				(euler_fluxes->A_1_2_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_2_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				(euler_fluxes->A_1_2_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_2_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );
				}

			}else{
						if(j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*e->ny[i]*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);			
						}else if(j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*e->ny[i];
						}else if(j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*e->ny[i];
						}else if(j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*e->ny[i];
						}
			}
		
    }

		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;
    return (-result);
 }

    double EulerBoundaryBilinearForm_vel_y::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerBoundaryBilinearForm_vel_y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerBoundaryBilinearForm_vel_y::clone() const { return new EulerBoundaryBilinearForm_vel_y(kappa, this->j,this->euler_fluxes, this->riemann_invariants, this->mirror_condition); }


//energy
  template<typename Real, typename Scalar>
    Scalar EulerBoundaryBilinearForm_e::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
 {
		int bdry =0;
		double* ghost_state = new double[4];
		double* dudu_j = new double[4];
		double* A_n = new double[4];
		bool solid = true;
		double constant =0.5;

    Scalar result = Scalar(0);	
    for (int i = 0;i < n;i++) 
    {				

			if((mirror_condition==true)||(solid==false)){ 
				bdry =riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


				if(bdry!=1){
							calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 3, A_n); //3te-Zeile A_n
					if(bdry!=2) 
							riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
				}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] *constant
        * euler_fluxes->A_1_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_1_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * euler_fluxes->A_2_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				(euler_fluxes->A_1_3_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_3_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				(euler_fluxes->A_1_3_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_3_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				(euler_fluxes->A_1_3_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_3_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				(euler_fluxes->A_1_3_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 euler_fluxes->A_2_3_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );

				}
				
     }
		}
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;
    return (-result);
  }

    double EulerBoundaryBilinearForm_e::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerBoundaryBilinearForm_e::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerBoundaryBilinearForm_e::clone() const { return new EulerBoundaryBilinearForm_e(kappa, this->j,this->euler_fluxes, this->riemann_invariants, this->mirror_condition); }



//--Linearforms---
//rho
 template<typename Real, typename Scalar>
    Scalar EulerBoundary_rho::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  		 
		 {
		double rho_new= 0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
		double lambda =0.;
		double* new_variables = new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
		int boundary; bool solid = false;

  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++){
			rho = ext[0]->val[i]; 		
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i];  
			rho_energy = ext[3]->val[i]; 


				 solid =true;  //solid wall 

     if((mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];

	//determine free-stream values
				riemann_invariants->get_free_stream(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], e->tx[i],e->ty[i],	new_variables, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext,boundary,solid);


					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];		

		      result += wt[i] *e->nx[i]*v->val[i]*0.5* (rho_new
		      * euler_fluxes->A_1_0_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
					+	ext[0]->val[i]
		      * euler_fluxes->A_1_0_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy)) ;
		      result += wt[i]* e->ny[i]*v->val[i]*0.5*( rho_new
		      * euler_fluxes->A_2_0_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new)  
						+ext[0]->val[i]
		      * euler_fluxes->A_2_0_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy));
		      result += wt[i]  * e->nx[i]*v->val[i] *0.5* (rho_v_x_new
		      * euler_fluxes->A_1_0_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
							+ ext[1]->val[i]
		      * euler_fluxes->A_1_0_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
		      result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_x_new
		      * euler_fluxes->A_2_0_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
		      + ext[1]->val[i]
		      * euler_fluxes->A_2_0_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
		      result += wt[i]* e->nx[i]*v->val[i] *0.5*(rho_v_y_new
		      * euler_fluxes->A_1_0_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
		      + ext[2]->val[i]
		      * euler_fluxes->A_1_0_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
		      result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_y_new
		      * euler_fluxes->A_2_0_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
		       +ext[2]->val[i]
		      * euler_fluxes->A_2_0_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
		      result += wt[i]* e->nx[i]*v->val[i] *0.5*( rho_energy_new
		      * euler_fluxes->A_1_0_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
		       +ext[3]->val[i]
		      * euler_fluxes->A_1_0_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy))  ;
		      result += wt[i]* e->ny[i]*v->val[i] * 0.5* ( rho_energy_new 
		      * euler_fluxes->A_2_0_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
		      + ext[3]->val[i] 
		      * euler_fluxes->A_2_0_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;

				if(boundary!=1)	
						result -= wt[i]*v->val[i]* 0.5 * 
								calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, 0);

	//lambda = calculate_lambda_max(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa);
//result -= wt[i]*v->val[i]* 0.5 *lambda*(rho_energy_new - ext[0]->val[i]);
				}
		
			
		}

		delete [] new_variables;

    return (-result);
    }      


    double EulerBoundary_rho::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    } 


    Ord EulerBoundary_rho::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* EulerBoundary_rho::clone() const
    {
					return new EulerBoundary_rho(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition);
    }

//rho_velocity_x
 template<typename Real, typename Scalar>
    Scalar EulerBoundary_v_x::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
		 {
		double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
		double lambda =0.;
		double* new_variables = new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
int boundary; bool solid = false;
  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++)       {
			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i]; 

				solid = true;

    if((mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];
			
				riemann_invariants->get_free_stream(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], e->tx[i],e->ty[i],	new_variables, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, boundary, solid);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];	

        result += wt[i] *e->nx[i]*v->val[i]*0.5* (rho_new
        * euler_fluxes->A_1_1_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
				+	ext[0]->val[i]
        * euler_fluxes->A_1_1_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy)) ;
        result += wt[i]* e->ny[i]*v->val[i]*0.5*( rho_new
        * euler_fluxes->A_2_1_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new)  
					+ext[0]->val[i]
        * euler_fluxes->A_2_1_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy));
        result += wt[i]  * e->nx[i]*v->val[i] *0.5* (rho_v_x_new
        * euler_fluxes->A_1_1_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
						+ ext[1]->val[i]
        * euler_fluxes->A_1_1_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_x_new
        * euler_fluxes->A_2_1_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[1]->val[i]
        * euler_fluxes->A_2_1_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*(rho_v_y_new
        * euler_fluxes->A_1_1_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[2]->val[i]
        * euler_fluxes->A_1_1_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_y_new
        * euler_fluxes->A_2_1_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[2]->val[i]
        * euler_fluxes->A_2_1_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*( rho_energy_new
        * euler_fluxes->A_1_1_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[3]->val[i]
        * euler_fluxes->A_1_1_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy))  ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* ( rho_energy_new 
        * euler_fluxes->A_2_1_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[3]->val[i] 
        * euler_fluxes->A_2_1_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;

			if(boundary!=1)		
						result -= wt[i]*v->val[i]* 0.5 * 
								calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, 1);


//lambda = calculate_lambda_max(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa);
//result -= wt[i]*v->val[i]* 0.5 *lambda*(rho_energy_new - ext[1]->val[i]);

				}else{//solid wall ->no mirror
							result += wt[i]*v->val[i]*e->nx[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
				}
			
			
		}
		delete [] new_variables;
     return (-result);
   }      


    double EulerBoundary_v_x::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    } 


    Ord EulerBoundary_v_x::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* EulerBoundary_v_x::clone() const
    {
					return new EulerBoundary_v_x(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition);
    }


//rho_velocity_y
 template<typename Real, typename Scalar>
    Scalar EulerBoundary_v_y::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
		 {
	double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
	double lambda =0.;
double* new_variables = new double[4];
double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
double rho, rho_v_x, rho_v_y, rho_energy;
int boundary; bool solid = false;

  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++)  {
			rho = ext[0]->val[i]; 
			rho_v_x = ext[1]->val[i];
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i]; 

					solid = true;

    if((mirror_condition==true)||(solid==false)){ 

			rho_ext = ext[4]->val[i];
			rho_v_x_ext = ext[5]->val[i];
			rho_v_y_ext = ext[6]->val[i];
			rho_energy_ext = ext[7]->val[i];

			riemann_invariants->get_free_stream(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], e->tx[i],e->ty[i],	new_variables, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, boundary, solid);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];			

        result += wt[i] *e->nx[i]*v->val[i]*0.5* (rho_new
        * euler_fluxes->A_1_2_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
				+	ext[0]->val[i]
        * euler_fluxes->A_1_2_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy)) ;
        result += wt[i]* e->ny[i]*v->val[i]*0.5*( rho_new
        * euler_fluxes->A_2_2_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new)  
					+ext[0]->val[i]
        * euler_fluxes->A_2_2_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy));
        result += wt[i]  * e->nx[i]*v->val[i] *0.5* (rho_v_x_new
        * euler_fluxes->A_1_2_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
						+ ext[1]->val[i]
        * euler_fluxes->A_1_2_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_x_new
        * euler_fluxes->A_2_2_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[1]->val[i]
        * euler_fluxes->A_2_2_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*(rho_v_y_new
        * euler_fluxes->A_1_2_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[2]->val[i]
        * euler_fluxes->A_1_2_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_y_new
        * euler_fluxes->A_2_2_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[2]->val[i]
        * euler_fluxes->A_2_2_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*( rho_energy_new
        * euler_fluxes->A_1_2_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[3]->val[i]
        * euler_fluxes->A_1_2_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy))  ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* ( rho_energy_new 
        * euler_fluxes->A_2_2_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[3]->val[i] 
        * euler_fluxes->A_2_2_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;

			if(boundary!=1)		
						result -= wt[i]*v->val[i]* 0.5 * 
								calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, 2);


//lambda = calculate_lambda_max(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa);
//result -= wt[i]*v->val[i]* 0.5 *lambda*(rho_energy_new - ext[2]->val[i]);

			}else{//solid wall
				result += wt[i]*v->val[i]*e->ny[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
			}
		 
		}
		delete [] new_variables;
    return (-result);
   }     


    double EulerBoundary_v_y::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    } 


    Ord EulerBoundary_v_y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* EulerBoundary_v_y::clone() const
    {
					return new EulerBoundary_v_y(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition);
    }


//rho_energy
 template<typename Real, typename Scalar>
    Scalar EulerBoundary_e::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
		 {
	double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
	double lambda =0.;
	double* new_variables = new double[4];
	double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
	double rho, rho_v_x, rho_v_y, rho_energy;
	int boundary; bool solid = false;

  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++) {
		rho = ext[0]->val[i]; 
		rho_v_x = ext[1]->val[i]; 
		rho_v_y = ext[2]->val[i]; 
		rho_energy = ext[3]->val[i]; 

				solid = true;

    if((mirror_condition==true)||(solid==false)){ 

		rho_ext = ext[4]->val[i];
		rho_v_x_ext = ext[5]->val[i];
		rho_v_y_ext = ext[6]->val[i];
		rho_energy_ext = ext[7]->val[i];

			riemann_invariants->get_free_stream(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], e->tx[i],e->ty[i],	new_variables, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, boundary, solid);

				rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];			

        result += wt[i] *e->nx[i]*v->val[i]*0.5* (rho_new
        * euler_fluxes->A_1_3_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
				+	ext[0]->val[i]
        * euler_fluxes->A_1_3_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy)) ;
        result += wt[i]* e->ny[i]*v->val[i]*0.5*( rho_new
        * euler_fluxes->A_2_3_0<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new)  
					+ext[0]->val[i]
        * euler_fluxes->A_2_3_0<Scalar>(rho, rho_v_x, rho_v_y, rho_energy));
        result += wt[i]  * e->nx[i]*v->val[i] *0.5* (rho_v_x_new
        * euler_fluxes->A_1_3_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
						+ ext[1]->val[i]
        * euler_fluxes->A_1_3_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_x_new
        * euler_fluxes->A_2_3_1<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[1]->val[i]
        * euler_fluxes->A_2_3_1<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*(rho_v_y_new
        * euler_fluxes->A_1_3_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[2]->val[i]
        * euler_fluxes->A_1_3_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* (rho_v_y_new
        * euler_fluxes->A_2_3_2<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[2]->val[i]
        * euler_fluxes->A_2_3_2<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;
        result += wt[i]* e->nx[i]*v->val[i] *0.5*( rho_energy_new
        * euler_fluxes->A_1_3_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
         +ext[3]->val[i]
        * euler_fluxes->A_1_3_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy))  ;
        result += wt[i]* e->ny[i]*v->val[i] * 0.5* ( rho_energy_new 
        * euler_fluxes->A_2_3_3<Scalar>(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new) 
        + ext[3]->val[i] 
        * euler_fluxes->A_2_3_3<Scalar>(rho, rho_v_x, rho_v_y, rho_energy) ) ;

			if(boundary!=1)		
						result -= wt[i]*v->val[i]* 0.5 * 
								calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, 3);

				//lambda = calculate_lambda_max(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i], rho_new,rho_v_x_new, rho_v_y_new, rho_energy_new, kappa);
				//result -= wt[i]*v->val[i]* 0.5 *lambda*(rho_energy_new - ext[3]->val[i]);
				}
			
		}
	
		delete [] new_variables;
            return (-result);
    }      


    double EulerBoundary_e::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    } 


    Ord EulerBoundary_e::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* EulerBoundary_e::clone() const
    {
					return new EulerBoundary_e(kappa,this->euler_fluxes, this->riemann_invariants, this->mirror_condition);
    }


