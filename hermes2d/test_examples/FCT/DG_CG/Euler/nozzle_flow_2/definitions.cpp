#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(int num_of_equations): WeakForm<double>(num_of_equations), num_of_equations(num_of_equations)
	{
		for(int k =0; k<num_of_equations;k++)		
			add_matrix_form(new DefaultMatrixFormVol<double>(k,k)); 
		
	}


	    WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
      const_cast<EulerEquationsWeakForm_Mass*>(this)->warned_nonOverride = false;
      return new EulerEquationsWeakForm_Mass(this->num_of_equations);
    }
    

//_----------------------------------------------------------------------------------------------------------

	
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
	EulerS::EulerS(double kappa,MeshSharedPtr mesh,NumericalFlux* num_flux,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, bool mirror_condition,int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)),kappa(kappa), riemann_invariants(new RiemannInvariants(kappa)), mirror_condition(mirror_condition),
    prev_density(prev_density), prev_density_vel_x(prev_density_vel_x), prev_density_vel_y(prev_density_vel_y), prev_energy(prev_energy), rho_ext(rho_ext), v1_ext (v1_ext), v2_ext(v2_ext), energy_ext (energy_ext), mesh(mesh) 
	{

		for(int k =0; k<4;k++)
		{
			for(int i = 0; i<4;i++)
				add_matrix_form_surf(new EulerS::EulerBoundaryBilinearForm(kappa,i,k,num_flux));
			
			add_vector_form_surf(new EulerS::EulerBoundaryLinearform(kappa, k,num_flux));
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
    wf = new EulerS(this->kappa,this->mesh,this->num_flux,this->rho_ext, this->v1_ext, this->v2_ext, this->energy_ext, this->prev_density, this->prev_density_vel_x, this->prev_density_vel_y, this->prev_energy, this->mirror_condition,4);

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

double steigung(double x, double y)
{
if(y>0)
 return (-0.25*0.5*PI*(Hermes::sin(PI*x/2.)));
else
 return (0.25*0.5*PI*(Hermes::sin(PI*x/2.)));
}

double normalized(double x, double y)
{

	return Hermes::sqrt(x*x+y*y);
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
 double w_L[4], w_R[4];

double nx, ny, tx, ty;
double nx_ghost, ny_ghost, tx_ghost, ty_ghost;
double pressure= 2./3.;

  for (int i = 0;i < n;i++) 
  {	


nx = e->nx[i];
ny = e->ny[i];
tx = e->tx[i];
ty = e->ty[i];
nx_ghost = nx;
ny_ghost = ny;
tx_ghost = tx;
ty_ghost = ty;

if((e->x[i]>0.)&&(e->x[i]<4.))
{
double m = steigung(e->x[i],e->y[i]);
double norm = normalized(1,m);

 if(e->y[i] <0)
{
	tx_ghost = 1./norm;
	ty_ghost = m/norm;
}else 
{
	tx_ghost = -1./norm;
	ty_ghost = -m/norm;
}
	nx_ghost = ty_ghost; ny_ghost = -tx_ghost;

}
	
				if((e->x[i]<8.)&&(e->x[i]>-2.))	{	solid = true; bdry =0;}
				else if(e->x[i]==8.){ //continue;
					bdry=4;solid==false; //sub out
				}else if(e->x[i]==-2.){ //continue;
					bdry =3.;solid==false;  //sub in
				}else throw Hermes::Exceptions::Exception("boundary");

		 if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 


(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state(bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx_ghost,ny_ghost,tx_ghost,ty_ghost, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state, solid,pressure);

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
//if(solid==false) result =0;
      return (-result);
  }  




    Ord EulerS::EulerBoundaryBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerS::EulerBoundaryBilinearForm::clone() const
    {
					return new EulerS::EulerBoundaryBilinearForm(kappa,this->entry_i, this->entry_j, this->num_flux);
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
  double pressure = 2./3.;
     double w_L[4];
double nx, ny, tx, ty, nx_ghost, ny_ghost, tx_ghost, ty_ghost;
  for (int i = 0;i < n;i++) 
{


nx = e->nx[i];
ny = e->ny[i];
tx = e->tx[i];
ty = e->ty[i];
nx_ghost = nx;
ny_ghost = ny;
tx_ghost = tx;
ty_ghost = ty;

if((e->x[i]>0.)&&(e->x[i]<4.))
{
double m = steigung(e->x[i],e->y[i]);
double norm = normalized(1,m);

 if(e->y[i] <0)
{
	tx_ghost = 1./norm;
	ty_ghost = m/norm;
}else 
{
	tx_ghost = -1./norm;
	ty_ghost = -m/norm;
}
	nx_ghost = ty_ghost; ny_ghost = -tx_ghost;
}


			if((e->x[i]<8.)&&(e->x[i]>-2.))	{	solid = true; bdry =0; }
				else if(e->x[i]==8.){ 
						bdry=4; solid==false; 				
				}else if(e->x[i]==-2.){ 
						bdry =3.; solid==false;					
				}else throw Hermes::Exceptions::Exception("boundary");

			rho = ext[0]->val[i];  
			rho_v_x = ext[1]->val[i]; 
			rho_v_y = ext[2]->val[i]; 
			rho_energy = ext[3]->val[i];

    if(((static_cast<EulerS*>(wf))->mirror_condition==true)||(solid==false)){ 
				rho_ext = ext[4]->val[i];
				rho_v_x_ext = ext[5]->val[i];
				rho_v_y_ext = ext[6]->val[i];
				rho_energy_ext = ext[7]->val[i];
			

(static_cast<EulerS*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx_ghost,ny_ghost,tx_ghost,ty_ghost, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], new_variables, solid,pressure);

					rho_new=new_variables[0]; rho_v_x_new=new_variables[1]; rho_v_y_new=new_variables[2]; rho_energy_new=new_variables[3];
							 w_L[0] = ext[0]->val[i];
        w_L[1] = ext[1]->val[i];
        w_L[2] = ext[2]->val[i];
        w_L[3] = ext[3]->val[i];

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
								
	//result += wt[i] * (static_cast<EulerS*>(wf))->num_flux->numerical_flux_i(this->i,w_L , new_variables, e->nx[i], e->ny[i]) * v->val[i];	

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
					return new EulerS::EulerBoundaryLinearform(this->kappa, this->entry_i, this->num_flux);
    }




