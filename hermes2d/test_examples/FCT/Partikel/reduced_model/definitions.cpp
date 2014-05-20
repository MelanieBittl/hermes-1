#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(int num_of_equations): WeakForm<double>(num_of_equations), num_of_equations(num_of_equations)
	{

		for(int k =0; k<num_of_equations;k++)			
			add_matrix_form(new DefaultMatrixFormVol<double>(k,k)); 	 

	}


	    WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
    EulerEquationsWeakForm_Mass* wf;
    wf = new EulerEquationsWeakForm_Mass(5);

    return wf;
    }

    


//_---------------Matrix K -------------------------------------------------------------------------------------------

	EulerK::EulerK(double gamma, MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(gamma)),gamma(gamma), riemann_invariants(new RiemannInvariants(gamma)),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g),
prev_density_p(prev_density_p)
	{
//gas phase
	for(int k =0; k<4;k++)
	{	for(int i = 0; i<4;i++)			
			add_matrix_form(new EulerK::EulerEquationsBilinearForm(i,k,gamma));	
	 add_vector_form(new EulerK::EulerEquationsLinearForm(k,gamma));
	}
//particle phase
	for(int k =4; k<5;k++)
	{	for(int i = 4; i<5;i++)			
			add_matrix_form(new EulerK::EulerEquationsBilinearForm(i,k,gamma,true));	
	 add_vector_form(new EulerK::EulerEquationsLinearForm(k,gamma,true));
	}
    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p));

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
this->prev_density_p,5);

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
			if(particle)
			{

				double v_x_g = ext[1]->val[i]/ext[0]->val[i]; 
				double v_y_g = ext[2]->val[i]/ext[0]->val[i]; 

				result += wt[i] * u->val[i] *(v_x_g * v->dx[i]+  v_y_g  * v->dy[i]);	
			}else{	

				result += wt[i] * u->val[i] *
				( (static_cast<EulerK*>(wf))->euler_fluxes->A_g(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i],0,entry_i,entry_j) 
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
		if(particle)
		{	 double v_x_g = ext[1]->val[i]/ext[0]->val[i]; 
			double v_y_g = ext[2]->val[i]/ext[0]->val[i]; 
			result += wt[i] * ext[4]->val[i]*( v->dx[i]*  v_x_g + v->dy[i]*v_y_g);	  

		}else{
			for(int k =0; k<4;k++)
			{
				  result += wt[i] * ext[k]->val[i]*(v->dx[i]*
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

	EulerSource::EulerSource(double particle_density,MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, int num_of_equations): WeakForm<double>(num_of_equations),particle_density(particle_density),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g),
prev_density_p(prev_density_p)

	{
//gas/particle phase
	for(int k =0; k<5;k++)
	{	for(int i = 0; i<5;i++)			
			add_matrix_form(new EulerSource::EulerSourceBilinearForm(i,k));	
	 add_vector_form(new EulerSource::EulerSourceLinearForm(k));
	}

    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p));

	};


	EulerSource ::~EulerSource ()
	{
	};
	
	WeakForm<double>* EulerSource::clone() const
    {
    EulerSource* wf;
    wf = new EulerSource(this->particle_density, this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,this->prev_density_p,5);

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

      double result = 0.;
/*if(entry_i!=4) return 0;

		double material_density = (static_cast<EulerSource*>(wf))->particle_density;  
if(entry_j==4)    
		for (int i = 0;i < n;i++)
		{
			if((e->y[i]>= -0.14)&&(e->y[i]<= 0.15)&&(e->x[i]>1.2)&&(e->x[i]<1.3)) 
			{
				// double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
				//if(mach>0.11){
						double rho = ext[4]->val[i];
		
							result += wt[i] *u->val* v->val[i];
				//}
			}

		}*/

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

double material_density = (static_cast<EulerSource*>(wf))->particle_density;
	

      double result = 0.;
if(entry_i==4){
	for (int i = 0;i < n;i++)
		{
			if((e->y[i]> -0.05)&&(e->y[i]< 0.05)&&((e->x[i]>1.06)&&(e->x[i]<=1.08))) 
			{
				// double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
				//if(mach>0.11){
						double rho = ext[4]->val[i];
		
							result += wt[i] * v->val[i]*(material_density*0.1);
				//}
			}

				      
		}	
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





//---------------------------Boundary only------------------------
	EulerBoundary::EulerBoundary(double gamma,MeshFunctionSharedPtr<double>  rho_ext_g, MeshFunctionSharedPtr<double>  v1_ext_g, MeshFunctionSharedPtr<double>  v2_ext_g, MeshFunctionSharedPtr<double>  energy_ext_g,MeshFunctionSharedPtr<double>  rho_ext_p,
MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(gamma)),gamma(gamma), riemann_invariants(new RiemannInvariants(gamma)),
    prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g), 
	rho_ext_g(rho_ext_g), v1_ext_g(v1_ext_g), v2_ext_g(v2_ext_g), energy_ext_g(energy_ext_g),
   prev_density_p(prev_density_p), rho_ext_p(rho_ext_p)
	{
//gas phase
		for(int k =0; k<4;k++)
		{
			for(int i = 0; i<4;i++)
				add_matrix_form_surf(new EulerBoundary::EulerBoundaryBilinearForm(gamma,i,k));
			
			add_vector_form_surf(new EulerBoundary::EulerBoundaryLinearform(gamma, k));
		}

//particle phase
		for(int k =4; k<5;k++)
		{
			for(int i = 4; i<5;i++)
				add_matrix_form_surf(new EulerBoundary::EulerBoundaryBilinearForm(gamma,i,k,true));
			
			add_vector_form_surf(new EulerBoundary::EulerBoundaryLinearform(gamma, k,true));
		}

    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g, rho_ext_g, v1_ext_g, v2_ext_g, energy_ext_g, prev_density_p, rho_ext_p));

	};


	EulerBoundary ::~EulerBoundary ()
	{

		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	WeakForm<double>* EulerBoundary::clone() const
    {
    EulerBoundary* wf;
    wf = new EulerBoundary(this->gamma,this->rho_ext_g, this->v1_ext_g, this->v2_ext_g, this->energy_ext_g,
		this->rho_ext_p, 
 this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,
this->prev_density_p,5);

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
 return -(Hermes::sin(x*PI)*PI)/(2.*1.28/0.72+2.); 
else
 return (Hermes::sin(x*PI)*PI)/(2.*1.28/0.72+2.);
}

double normalized(double x, double y)
{

	return Hermes::sqrt(x*x+y*y);
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
double nx_ghost, ny_ghost, tx_ghost, ty_ghost;
	double nx,ny,tx, ty;
	double constant = 1.;
  	double result = 0.;
  for (int i = 0;i < n;i++) 
  {	

			nx = e->nx[i];
			ny = e->ny[i];
			tx = e->tx[i];
			ty = e->ty[i];
if(e->x[i]== 1.)
{
nx = -1.; ny = 0.;
ty = nx;
tx = -ny;

}
nx_ghost = nx;
ny_ghost = ny;
tx_ghost = tx;
ty_ghost = ty;

/*
if((e->x[i]<1.))
{
		double m = steigung(e->x[i],e->y[i]);
		double norm = normalized(1,m);

		 if(e->y[i] >0)
		{
			tx_ghost = -1./norm;
			ty_ghost = -m/norm;
		}else 
		{
			tx_ghost = 1./norm;
			ty_ghost = m/norm;
		}
	nx_ghost = ty_ghost; ny_ghost = -tx_ghost;

}*/

		bdry = 0;
		if(e->x[i]== 0.)
		{ 	//bdry = 4;
			//double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
			//if(mach>0.99) 
				bdry = 2;
		}
		if(e->x[i]== 2){ 
			//bdry =3;
			//double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
			//if(mach>0.99) 
				bdry = 1;
		}


		if(bdry==1) constant =1.;
		else constant = 0.5;

	//----------particle---------------------
	if(entry_i>=4)
		{	if(bdry==0) continue;

			rho = 	  ext[8]->val[i];    
			rho_ext = ext[9]->val[i];

			double v_x_g = ext[1]->val[i]/ext[0]->val[i]; 
			double v_y_g = ext[2]->val[i]/ext[0]->val[i]; 

			if((bdry ==3)||(bdry ==2)) //inlet
			{
				return 0.;
			}			

		result += wt[i] * u->val[i] *v->val[i]*( v_x_g * nx + v_y_g * ny);



	}else{//-------------gas phase-------------------
			rho = 		ext[0]->val[i];  
			rho_v_x =   ext[1]->val[i]; 
			rho_v_y = 	ext[2]->val[i]; 
			rho_energy= ext[3]->val[i];

    
				rho_ext = 		ext[4]->val[i];
				rho_v_x_ext = 	ext[5]->val[i];
				rho_v_y_ext = 	ext[6]->val[i];
				rho_energy_ext =ext[7]->val[i];

	 		if(bdry!=0){ 
				(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state(bdry,rho, rho_v_x, rho_v_y, rho_energy, nx_ghost,ny_ghost,tx_ghost,ty_ghost, rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext, ghost_state);

				if(bdry!=1){
						Boundary_helpers::calculate_A_n(rho, rho_v_x, rho_v_y, rho_energy, nx,ny , ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3], gamma, entry_i,A_n); //ite-Zeile A_n

					if(bdry!=2) 
					{
							(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_du_du(rho, rho_v_x, rho_v_y, rho_energy, nx,ny ,tx, ty, ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],bdry,entry_j,dudu_j);//j-te Spalte
					}
				}

				
		    result += wt[i] * u->val[i] *v->val[i]* constant*
		    ( (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,entry_j) 
		      * nx +
		     (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,entry_j)  
		      * ny);
						 	

				if (bdry==2)
			  		result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
				else if(bdry!=1){ 
				 	result += wt[i]*v->val[i] *u->val[i]*0.5* A_n[entry_j];
					for(int k =0;k<4;k++)
					{
						result += wt[i]*v->val[i] *u->val[i]*0.5*(
						((static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],0,entry_i,k) 
										  * nx+
						(static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(ghost_state[0], ghost_state[1], ghost_state[2],ghost_state[3],1,entry_i,k) 
										  * ny -A_n[k])* dudu_j[k]);
					} 
				}
			
			}else{ 
					if(entry_i==1)
					{
							if(entry_j==0){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*0.5*nx*
							(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);	
	
							}else if(entry_j==1){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[1]->val[i]/ext[0]->val[i]*nx;
							}else if(entry_j==2){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[2]->val[i]/ext[0]->val[i]*nx;
							}else if(entry_j==3){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*nx;
							}
					}else if(entry_i==2)
					{
							if(entry_j==0){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*0.5*ny*
							(ext[1]->val[i]*ext[1]->val[i]+ext[2]->val[i]*ext[2]->val[i])/(ext[0]->val[i]*ext[0]->val[i]);			
							}else if(entry_j==1){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[1]->val[i]/ext[0]->val[i]*ny;
							}else if(entry_j==2){
								result -= wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ext[2]->val[i]/ext[0]->val[i]*ny;
							}else if(entry_j==3){
								result += wt[i]*v->val[i] *u->val[i]*(gamma-1.)*ny;
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
double nx,ny,tx, ty;
double nx_ghost, ny_ghost, tx_ghost, ty_ghost;
int bdry; 
  double result = 0.;
  for (int i = 0;i < n;i++) 
	{
		bdry = 0;
		if(e->x[i]== 0.)
		{ 	//bdry = 4;
			//double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
			//if(mach>0.99) 
				bdry = 2;
		}
		if(e->x[i]== 2){ 
			//bdry =3;
			//double mach = QuantityCalculator::calc_mach(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i], 1.4);
			//if(mach>0.99) 
				bdry = 1;
		}


			nx = e->nx[i];
			ny = e->ny[i];
			tx = e->tx[i];
			ty = e->ty[i];
if(e->x[i]== 1.)
{
nx = -1.; ny = 0.;
ty = nx;
tx = -ny;

}
nx_ghost = nx;
ny_ghost = ny;
tx_ghost = tx;
ty_ghost = ty;


/*
if((e->x[i]<1.))
{
double m = steigung(e->x[i],e->y[i]);
double norm = normalized(1,m);

 if(e->y[i] >0)
{
	tx_ghost = -1./norm;
	ty_ghost = -m/norm;
}else 
{
	tx_ghost = 1./norm;
	ty_ghost = m/norm;
}
	nx_ghost = ty_ghost; ny_ghost = -tx_ghost;

}*/



	//----------particle---------------------
	if(entry_i==4)
	{	
		if(bdry==0) continue;	

			rho = ext[8]->val[i];    
			rho_ext = ext[9]->val[i];

			double v_x_ext = ext[5]->val[i]/ext[4]->val[i]; 
			double v_y_ext = ext[6]->val[i]/ext[4]->val[i]; 

			double v_x_g = ext[1]->val[i]/ext[0]->val[i]; 
			double v_y_g = ext[2]->val[i]/ext[0]->val[i]; 

			if((bdry ==3)||(bdry ==2)) //inlet
			{
					result += wt[i] *v->val[i]*( v_x_ext * nx + v_y_ext * ny)*rho_ext;	
			}else 
				result += wt[i] *v->val[i]*( v_x_g * nx + v_y_g * ny)*rho;			

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
			(static_cast<EulerBoundary*>(wf))->riemann_invariants->get_ghost_state( bdry,ext[0]->val[i], ext[1]->val[i], ext[2]->val[i],ext[3]->val[i], nx_ghost,ny_ghost,tx_ghost,ty_ghost, ext[4]->val[i], ext[5]->val[i], ext[6]->val[i],ext[7]->val[i], ghost_state);

		rho_new=ghost_state[0]; 
		rho_v_x_new=ghost_state[1]; 
		rho_v_y_new=ghost_state[2]; 
		rho_energy_new=ghost_state[3];	
				for(int k =0; k<4;k++)
				{
				  result += wt[i] *nx*v->val[i]*0.5* (ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,0,entry_i,k) 
						+	ext[k]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,0,entry_i,k)) ;
				  result += wt[i]* ny*v->val[i]*0.5*( ghost_state[k]*
				   (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new,1,entry_i,k) 
							+ext[k]->val[i]
				  * (static_cast<EulerBoundary*>(wf))->euler_fluxes->A_g(rho, rho_v_x, rho_v_y, rho_energy,1,entry_i,k)) ;

				}

				if(bdry!=1)		
					result -= wt[i]*v->val[i]* 0.5 * 
										Boundary_helpers::calculate_A_n_U(rho, rho_v_x, rho_v_y, rho_energy, nx, ny,  rho_new, rho_v_x_new, rho_v_y_new, rho_energy_new, gamma, entry_i);
		}else{//solid wall
				if(entry_i==1)
					result += wt[i]*v->val[i]*nx*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
				else if(entry_i==2)
					result += wt[i]*v->val[i]*ny*QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y, rho_energy, gamma);
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


//--------------------------Penalty Term------------------------
	EulerPenalty::EulerPenalty(double sigma,double particle_density,MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, double eps,int num_of_equations): WeakForm<double>(num_of_equations), sigma(sigma), eps(eps), prev_density_g(prev_density_g), prev_density_vel_x_g(prev_density_vel_x_g), prev_density_vel_y_g(prev_density_vel_y_g), prev_energy_g(prev_energy_g), particle_density(particle_density),
prev_density_p(prev_density_p)
	{
//gas phase
		for(int k =0; k<4;k++)
		{
			for(int i = 0; i<4;i++)
				add_matrix_form_surf(new EulerPenalty::PenaltyBilinearForm(sigma,eps,i,k));
			
			add_vector_form_surf(new EulerPenalty::PenaltyLinearForm(sigma, k));
		}

//particle phase
		for(int k =4; k<5;k++)
		{
			for(int i = 4; i<5;i++)
				add_matrix_form_surf(new EulerPenalty::PenaltyBilinearForm(sigma,eps,i,k,true));
			
			add_vector_form_surf(new EulerPenalty::PenaltyLinearForm(sigma, k,true));
		}

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p));

	};


	
	WeakForm<double>* EulerPenalty::clone() const
    {
    EulerPenalty* wf;
    wf = new EulerPenalty(this->sigma,this->particle_density,this->prev_density_g, this->prev_density_vel_x_g, this->prev_density_vel_y_g, this->prev_energy_g,
this->prev_density_p,this->eps,8);

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


//-------Penalty Bilinearforms---------
   double EulerPenalty::PenaltyBilinearForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
{
  double result = 0.;
  for (int i = 0;i < n;i++) 
  {		
		
if((e->x[i]== 0) ||(e->x[i]== 1)) continue;
if((entry_i==0) ||(entry_i==3)||(entry_i==4)) return 0.; 
if((entry_j==0) ||(entry_j==3)||(entry_j==4) ) return 0.;
  
		if(particle){


		}else{
		double vn = ext[1]->val[i]/ext[0]->val[i]*e->nx[i] + ext[2]->val[i]/ext[0]->val[i]*e->ny[i];
		double factor = ((2.*vn*vn+eps)/(std::sqrt(vn*vn+eps)));


				if(entry_i==1)
				{
						if(entry_j==1){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(e->nx[i]);
						}else if(entry_j==2){
							result += wt[i]*v->val[i] *u->val[i]*e->ny[i]*e->nx[i];
						}
				}else if(entry_i==2)
				{
					 if(entry_j==1){
							result += wt[i]*v->val[i] *u->val[i]*e->ny[i]*e->nx[i];
						}else if(entry_j==2){
							result += wt[i]*v->val[i] *u->val[i]*Hermes::sqr(e->ny[i]);
						}
				}
		}
			
    }


      return sigma*(-result);
  }  




    Ord EulerPenalty::PenaltyBilinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    MatrixFormSurf<double>* EulerPenalty::PenaltyBilinearForm::clone() const
    {
					return new EulerPenalty::PenaltyBilinearForm(sigma,eps,this->entry_i, this->entry_j, this->particle);
    }





//------------------------------------------------------------------
//----------------------------Linearform Penalty---------------------------------
//----------------------------------------------


    double EulerPenalty::PenaltyLinearForm::value(int n, double *wt, Func<double> *u_ext[],  Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
 {
double material_density = (static_cast<EulerPenalty*>(wf))->particle_density;
  double result = 0.;
  for (int i = 0;i < n;i++) 
	{
if((e->x[i]== 0) ||(e->x[i]== 1)) continue; 
if((entry_i==0) ||(entry_i==3)||(entry_i==4)) return 0.; 
		double alpha_p  =ext[4]->val[i]/material_density ;
			double alpha_g = 1.-alpha_p;

		if(particle)
		{


		}else{
			double vn = (ext[1]->val[i]*e->nx[i] + ext[2]->val[i]*e->ny[i])/ext[0]->val[i];
			double rho =  ext[0]->val[i]/alpha_g;
				if(entry_i==1)
					result += wt[i]*v->val[i]*e->nx[i]*vn*std::fabs(vn)*rho*rho;
				else if(entry_i==2)
					result += wt[i]*v->val[i]*e->ny[i]*vn*std::fabs(vn)*rho*rho;
		}			
	 }

     return sigma*(-result);
  }      



    Ord EulerPenalty::PenaltyLinearForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(10);
    }

    VectorFormSurf<double>* EulerPenalty::PenaltyLinearForm::clone() const
    {
					return new EulerPenalty::PenaltyLinearForm(this->sigma, this->entry_i, this->particle);
    }














