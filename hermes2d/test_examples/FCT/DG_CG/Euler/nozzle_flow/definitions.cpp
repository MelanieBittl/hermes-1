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

//_----------------------------------------------------------------------------------------------------------

	EulerKS::EulerKS(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), rho_ext(rho_ext), v1_ext (v1_ext), v2_ext(v2_ext), energy_ext (energy_ext) 
	{

for(int k =0; k<4;k++)
{
		for(int i = 0; i<4;i++)
		{	add_matrix_form_surf(new EulerKS::EulerBoundaryBilinearForm(kappa,i,k));
			add_matrix_form(new EulerKS::EulerEquationsBilinearForm(i,k,kappa));

		}
		add_vector_form_surf(new EulerKS::EulerBoundaryLinearform(kappa, k));
}

    
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


//---------------------K----------------

    double EulerKS::EulerEquationsBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		 {
      double result = 0.;
      for (int i = 0;i < n;i++) 
      {

        result += wt[i] * u->val[i] *
        ( (static_cast<EulerKS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
          * v->dx[i]+        
         (static_cast<EulerKS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,entry_j) 
          * v->dy[i]);
		
      }
      return result;
    }
    


    Ord EulerKS::EulerEquationsBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormVol<double>* EulerKS::EulerEquationsBilinearForm::clone() const
    {
    return new EulerKS::EulerEquationsBilinearForm(this->entry_i, this->entry_j, this->kappa);
    }


//---Boudary Bilinearforms---------
   double EulerKS::EulerBoundaryBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
	int bdry =0;
	double* ghost_state = new double[4];for(int l=0;l<4;l++) ghost_state[l]=0.;
	double * A_n = new double[4];
	double* dudu_j= new double[4];
	bool solid = false;
	double constant = 1.;
  double result = 0.;
  for (int i = 0;i < n;i++) 
//if(e->x[i]!=-2.)
if((e->x[i]<8.)&&(e->x[i]>-2.))
  {		
if((e->x[i]<8.)&&(e->x[i]>-2.))	{	solid = true; bdry =0;}
else solid = false;

		 if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 

				if((e->x[i]<8.)&&(e->x[i]>-2.))	{	solid = true; bdry =0;}
				else if(e->x[i]==8.){ bdry=3;}
				else if(e->x[i]==-2.) bdry =4.;
				else throw Hermes::Exceptions::Exception("boundary");

(static_cast<EulerKS*>(wf))->riemann_invariants->get_ghost_state(bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid);

if(bdry==1) constant =1.;
else constant = 0.5;


			if(bdry!=1){
					Boundary_helpers::calculate_A_n(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], kappa, entry_i,A_n); //ite-Zeile A_n

				if(bdry!=2) 
				{
						(static_cast<EulerKS*>(wf))->riemann_invariants->get_du_du(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i] ,e->tx[i], e->ty[i], ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,entry_j,dudu_j);//j-te Spalte
				}
			}

				
        result += wt[i] * u->val[i] *v->val[i]* constant*
        ( (static_cast<EulerKS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
          * e->nx[i] +
         (static_cast<EulerKS*>(wf))->euler_fluxes->A(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],1,entry_i,entry_j)  
          * e->ny[i]);
			

			if (bdry==2)
	      		result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[entry_j];
			else if(bdry!=1){ 
			 	result += wt[i]*v->val[i] *u->val[i]*0.5*A_n[entry_j];
				for(int k =0;k<4;k++)
				{result += wt[i]*v->val[i] *u->val[i]*0.5*(
					((static_cast<EulerKS*>(wf))->euler_fluxes->A(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],0,entry_i,k) 
									  * e->nx[i]+
					(static_cast<EulerKS*>(wf))->euler_fluxes->A(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],1,entry_i,k) 
									  * e->ny[i]  -A_n[k])* dudu_j[k]);
				} 
			}
			
		 }else{
				if(entry_i==1)
				{
						if(entry_j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*e->nx[i]*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);		
	
						}else if(entry_j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*e->nx[i];
						}else if(entry_j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*e->nx[i];
						}else if(entry_j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*e->nx[i];
						}
				}else if(entry_i==2)
				{
						if(entry_j==0){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*0.5*e->ny[i]*
						(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);			
						}else if(entry_j==1){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[1]->val[i]/ext[0]->val[i]*e->ny[i];
						}else if(entry_j==2){
							result -= wt[i]*v->val[i] *u->val[i]*(kappa-1.)*ext[2]->val[i]/ext[0]->val[i]*e->ny[i];
						}else if(entry_j==3){
							result += wt[i]*v->val[i] *u->val[i]*(kappa-1.)*e->ny[i];
						}

				}
			}
    }
		delete [] ghost_state;
		delete [] A_n;
		delete [] dudu_j;

      return (-result);
  }  




    Ord EulerKS::EulerBoundaryBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerKS::EulerBoundaryBilinearForm::clone() const
    {
					return new EulerKS::EulerBoundaryBilinearForm(kappa,this->entry_i, this->entry_j);
    }





//------------------------------------------------------------------
//----------------------------Linearforms---------------------------------
//----------------------------------------------


    double EulerKS::EulerBoundaryLinearform::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
		 {
		double rho_new=0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
		double lambda =0.;
		double* new_variables = new double[4];
		double rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext;
		double rho, rho_v_x, rho_v_y, rho_energy;
int bdry; bool solid = false;
  double result = 0.;
  for (int i = 0;i < n;i++) 
//if(e->x[i]==-2.)
if((e->x[i]==8.)||(e->x[i]==-2.))
{

if((e->x[i]<8.)&&(e->x[i]>-2.))	{	solid = true; bdry =0;}
else solid = false;
if(e->x[i]==-2.) bdry =3;
else bdry = 4.;

			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];

    if(((static_cast<EulerKS*>(wf))->mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];
			

(static_cast<EulerKS*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], e->nx[i],e->ny[i],e->tx[i],e->ty[i], ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], new_variables, solid);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];	

				for(int k =0; k<4;k++)
				{
							  result += wt[i] *e->nx[i]*v->val[i]*0.5* (new_variables[k]*
							   (static_cast<EulerKS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,entry_i,k) 
									+	ext[k]->val[i]
							  * (static_cast<EulerKS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,k)) ;
							  result += wt[i]* e->ny[i]*v->val[i]*0.5*( new_variables[k]*
							   (static_cast<EulerKS*>(wf))->euler_fluxes->A(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,entry_i,k) 
										+ext[k]->val[i]
							  * (static_cast<EulerKS*>(wf))->euler_fluxes->A(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,k)) ;

				}

		if(bdry!=1)		
			result -= wt[i]*v->val[i]* 0.5 * 
								Boundary_helpers::calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, e->nx[i], e->ny[i],  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, kappa, entry_i);


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



    Ord EulerKS::EulerBoundaryLinearform::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* EulerKS::EulerBoundaryLinearform::clone() const
    {
					return new EulerKS::EulerBoundaryLinearform(this->kappa, this->entry_i);
    }



