#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(int num_of_equations): WeakForm<double>(num_of_equations), num_of_equations(num_of_equations)
	{
    add_matrix_form(new EulerEquationsBilinearFormTime(0));  //density
    add_matrix_form(new EulerEquationsBilinearFormTime(1));		//density_vel_x
    add_matrix_form(new EulerEquationsBilinearFormTime(2));		//density_vel_y
    add_matrix_form(new EulerEquationsBilinearFormTime(3));		//energy
	}


	    WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
      const_cast<EulerEquationsWeakForm_Mass*>(this)->warned_nonOverride = false;
      return new EulerEquationsWeakForm_Mass(this->num_of_equations);
    }
    
    
    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormTime::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double EulerEquationsBilinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormTime::clone() const { return new EulerEquationsBilinearFormTime(component_i); }



	EulerKS::EulerKS(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), rho_ext(rho_ext), v1_ext (v1_ext), v2_ext(v2_ext), energy_ext (energy_ext) 
	{

 		add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_rho(kappa,0));
 		add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_rho(kappa,1));
 		add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_rho(kappa,2));
 		add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_rho(kappa,3));

    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_x( kappa,0));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_x( kappa,1));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_x( kappa,2));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_x( kappa,3));

    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_y( kappa,0));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_y( kappa,1));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_y( kappa,2));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_vel_y( kappa,3));

    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_e( kappa,0));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_e( kappa,1));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_e( kappa,2));
    add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm_e( kappa,3));

    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensity(0));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensity(1));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensity(2));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensity(3));

   add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelX(kappa,0));
   add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelX(kappa,1));
   add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelX(kappa,2));
   add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelX(kappa,3));

    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelY(kappa,0));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelY(kappa,1));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelY(kappa,2));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormDensityVelY(kappa,3));

    add_matrix_form(new EulerKS::EulerEquationsBilinearFormEnergy(kappa,0));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormEnergy(kappa,1));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormEnergy(kappa,2));
    add_matrix_form(new EulerKS::EulerEquationsBilinearFormEnergy(kappa,3));
    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, rho_ext, v1_ext, v2_ext, energy_ext));

	};


	EulerKS ::~EulerKS ()
	{

		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerKS::clone() const
    {
    EulerKS* wf;
    wf = new EulerKS(this->kappa,this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->mirror_condition,4);

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



  template<typename Real, typename Scalar>
    Scalar EulerKS::EulerBoundaryBilinearForm_rho::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
{
	int bdry =0;
	double* ghost_state = new double[4];
	double * A_n = new double[4];
	double* dudu_j= new double[4];
	bool solid = true;
	double constant = 0.5;
  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++) 
  {		


		 if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 
				bdry =(static_cast<EulerKS*>(wf))->riemann_invariants->get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);

			if(bdry!=1){
						Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 0, A_n); //0te-Zeile A_n
				if(bdry!=2) 
						(static_cast<EulerKS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
			}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}

				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );

				}
			
		 }
    }
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;

      return (-result);
  }  


    double EulerKS::EulerBoundaryBilinearForm_rho::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerKS::EulerBoundaryBilinearForm_rho::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerKS::EulerBoundaryBilinearForm_rho::clone() const
    {
					return new EulerKS::EulerBoundaryBilinearForm_rho(kappa, this->j);
    }

//vel_x
    template<typename Real, typename Scalar>
    Scalar EulerKS::EulerBoundaryBilinearForm_vel_x::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
{
	int bdry =0;
	double* ghost_state = new double[4];
	double * A_n = new double[4];
	double* dudu_j= new double[4];
	bool solid = true;
	double constant = 0.5;
  Scalar result = Scalar(0);
  for (int i = 0;i < n;i++) 
  {			


			if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 

	bdry =(static_cast<EulerKS*>(wf))->riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


	if(bdry!=1){
				Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 1, A_n); //1te-Zeile A_n
		if(bdry!=2) 
				(static_cast<EulerKS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
	}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] *constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
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

    double EulerKS::EulerBoundaryBilinearForm_vel_x::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerKS::EulerBoundaryBilinearForm_vel_x::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerKS::EulerBoundaryBilinearForm_vel_x::clone() const { return new EulerKS::EulerBoundaryBilinearForm_vel_x(kappa, this->j); }

//vel_y
  template<typename Real, typename Scalar>
    Scalar EulerKS::EulerBoundaryBilinearForm_vel_y::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
		int bdry =0;
		double* ghost_state = new double[4];
		double* dudu_j = new double[4];
		double* A_n = new double[4];
	bool solid = true;
	double constant = 0.5;

    Scalar result = Scalar(0);
    for (int i = 0;i < n;i++) 
    {				


			if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 
					bdry =(static_cast<EulerKS*>(wf))->riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


					if(bdry!=1){
								Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 2, A_n); //2te-Zeile A_n
						if(bdry!=2) 
								(static_cast<EulerKS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
					}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				 if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
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

    double EulerKS::EulerBoundaryBilinearForm_vel_y::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerKS::EulerBoundaryBilinearForm_vel_y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerKS::EulerBoundaryBilinearForm_vel_y::clone() const { return new EulerKS::EulerBoundaryBilinearForm_vel_y(kappa, this->j); }


//energy
  template<typename Real, typename Scalar>
    Scalar EulerKS::EulerBoundaryBilinearForm_e::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
 {
		int bdry =0;
		double* ghost_state = new double[4];
		double* dudu_j = new double[4];
		double* A_n = new double[4];
	bool solid = true;
	double constant = 0.5;

    Scalar result = Scalar(0);	
    for (int i = 0;i < n;i++) 
    {				

			if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 
				bdry =(static_cast<EulerKS*>(wf))->riemann_invariants->	get_bdry_info(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);


				if(bdry!=1){
							Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, 3, A_n); //3te-Zeile A_n
					if(bdry!=2) 
							(static_cast<EulerKS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,j,dudu_j);//j-te Spalte
				}

				if(j==0){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->ny[i]*v->val[i];
				}else if(j==1){
        result += wt[i] * u->val[i] *constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}else if(j==2){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * e->ny[i]*v->val[i];
				}else if(j==3){
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] * constant
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
				}
				if (bdry==2)
		      result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				else if(bdry!=1){ 
				 result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[j];
				result += wt[i]*v->val[i] *u->val[i]*0.5*(
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]  -A_n[0])* dudu_j[0] +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[1])*dudu_j[1]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[2])*dudu_j[2]  +
				((static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->nx[i]+
				 (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3]) 
								  * e->ny[i]-A_n[3])*dudu_j[3] );

				}
				
     }
		}
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;
    return (-result);
  }

    double EulerKS::EulerBoundaryBilinearForm_e::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerKS::EulerBoundaryBilinearForm_e::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerKS::EulerBoundaryBilinearForm_e::clone() const { return new EulerKS::EulerBoundaryBilinearForm_e(kappa, this->j); }



//---------------------K----------------

    template<typename Real, typename Scalar>
    Scalar EulerKS::EulerEquationsBilinearFormDensity::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
		 {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==1){
				
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }
    
     

    double EulerKS::EulerEquationsBilinearFormDensity::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerKS::EulerEquationsBilinearFormDensity::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerKS::EulerEquationsBilinearFormDensity::clone() const
    {
    return new EulerKS::EulerEquationsBilinearFormDensity(this->j);
    }



    template<typename Real, typename Scalar>
    Scalar EulerKS::EulerEquationsBilinearFormDensityVelX::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerKS::EulerEquationsBilinearFormDensityVelX::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerKS::EulerEquationsBilinearFormDensityVelX::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerKS::EulerEquationsBilinearFormDensityVelX::clone() const { return new EulerKS::EulerEquationsBilinearFormDensityVelX(kappa,this->j); }



    template<typename Real, typename Scalar>
    Scalar EulerKS::EulerEquationsBilinearFormDensityVelY::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerKS::EulerEquationsBilinearFormDensityVelY::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerKS::EulerEquationsBilinearFormDensityVelY::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerKS::EulerEquationsBilinearFormDensityVelY::clone() const { return new EulerKS::EulerEquationsBilinearFormDensityVelY(kappa,this->j); }


  template<typename Real, typename Scalar>
    Scalar EulerKS::EulerEquationsBilinearFormEnergy::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			 if (j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerKS*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerKS::EulerEquationsBilinearFormEnergy::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerKS::EulerEquationsBilinearFormEnergy::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(20);
    }


    MatrixFormVol<double>* EulerKS::EulerEquationsBilinearFormEnergy::clone() const { return new EulerKS::EulerEquationsBilinearFormEnergy(kappa,this->j); }






