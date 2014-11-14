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
	EulerS::EulerS(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,  NumericalFlux* num_flux, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),num_flux(num_flux),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), rho_ext(rho_ext), v1_ext (v1_ext), v2_ext(v2_ext), energy_ext (energy_ext) 
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
    wf = new EulerS(this->kappa,this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy,this->num_flux, this->mirror_condition,4);

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
double nx, ny, tx, ty;
  for (int i = 0;i < n;i++) 
  {		
		if((e->x[i]>= (1./6.))&&(e->y[i]==0)){
	 			solid =true;  //solid wall 
	 			bdry=0.;
	 	}
nx = e->nx[i];
ny = e->ny[i];
tx = e->tx[i];
ty = e->ty[i];

		 if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 
				if((e->x[i]>= (1./6.))&&(e->y[i]==0)){
					solid =true;  //solid wall 
					bdry=0.;
				}else if((e->x[i]< (1./6.))&&(e->y[i]==0)) bdry=1;
				else if((e->x[i]==4)&&(e->y[i]>0)) bdry=1;	
				else {
			//bdry = (static_cast<EulerS*>(wf))->riemann_invariants->get_bdry_info_short(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], solid);		

				//	if((bdry==5)||(bdry==3)) 
						//bdry=2;
					//if(bdry==4) bdry =1;
				
				if(ext[5]->val[i]==0) bdry =1;
				else bdry=2;
				}


(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state(bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx,ny,tx,ty, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);

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
if(solid==false) result = 0;
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
//----------------------------Linearforms---------------------------------
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
   double w_L[4];
  for (int i = 0;i < n;i++) 
{
	if((e->x[i]>= (1./6.))&&(e->y[i]==0)){
	 			solid =true;  //solid wall 
	 			bdry=0.;
	 	}

			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];

    if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];
				if((e->x[i]>= (1./6.))&&(e->y[i]==0)){
					solid =true;  //solid wall 
					bdry=0.;
				}else if((e->x[i]< (1./6.))&&(e->y[i]==0)) bdry=1;
				else if((e->x[i]==4)&&(e->y[i]>0)) bdry=1;	
				else {
	/*bdry = (static_cast<EulerS*>(wf))->riemann_invariants->get_bdry_info_short(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], solid);		
			
					if((bdry==5)||(bdry==3)) 
						bdry=2;
					if(bdry==4) bdry =1;*/
	
					if(ext[5]->val[i]==0) bdry =1;
				else bdry=2;
				}

		

(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], new_variables, solid);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];
		 w_L[0] = ext[0]->val[i];
        w_L[1] = ext[1]->val[i];
        w_L[2] = ext[2]->val[i];
        w_L[3] = ext[3]->val[i];

				for(int k =0; k<4;k++)
				{
							  result += wt[i] *e->nx[i]*v->val[i]*0.5* (new_variables[k]*
							   (static_cast<EulerS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,entry_i,k) 
									+	ext[k]->val[i]
							  * (static_cast<EulerS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,k)) ;
							  result += wt[i]* e->ny[i]*v->val[i]*0.5*( new_variables[k]*
							   (static_cast<EulerS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,entry_i,k) 
										+ext[k]->val[i]
							  * (static_cast<EulerS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,k)) ;

				}

		if(bdry!=1)		
			result -= wt[i]*v->val[i]* 0.5 * 
								Boundary_helpers::calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, entry_i);

			
	//	result += wt[i] * (static_cast<EulerS*>(wf))->num_flux->numerical_flux_i(this->i,w_L , new_variables, e->nx[i], e->ny[i]) * v->val[i];	

				}else{//solid wall ->no mirror
						if(entry_i==1)
							result += wt[i]*v->val[i]*e->nx[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
						else if(entry_i==2)
							result += wt[i]*v->val[i]*e->ny[i]*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, kappa);
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

