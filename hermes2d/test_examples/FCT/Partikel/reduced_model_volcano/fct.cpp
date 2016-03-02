
//Assemble antidiffusive fluxes & Limiter
void antidiffusiveFlux(CSCMatrix<double>* mass_matrix,CSCMatrix<double>* lumped_matrix,CSCMatrix<double>* diffusion,CSCMatrix<double>* matrix_L,
SimpleVector<double>*  vec_bdry,SimpleVector<double>*  vec_source,double* u_L, double* flux_scalar, double* P_plus, double* P_minus, double* Q_plus, double* Q_minus, double* R_plus, double* R_minus, int dof_rho, double time_step, double gamma )
{ 
 int dof_v_x = dof_rho; int dof_v_y = dof_rho;  int dof_e = dof_rho; 
int ndof_g = dof_rho*4.;
	int ndof = diffusion->get_size();
	double* dt_u_L = new double[ndof];
	double alpha,f, plus, minus,mass, diff,lumped;
int dof_rhoE_start = dof_rho + dof_v_x + dof_v_y;
int dof_v_y_start = dof_rho + dof_v_x ;
	double f_rho_ij, f_rho_ji, alpha_rho;
	double f_p_ij, f_p_ji, alpha_p;
	double f_v_x_ij, f_v_x_ji, alpha_v_x;
	double f_v_y_ij, f_v_y_ji, alpha_v_y;

	double f_rho, f_rho_v_x, f_rho_v_y, f_rhoE;
	double f_v_x,f_v_y,f_p, v_x_i, v_x_j, v_y_i, v_y_j, p_i, p_j;

	double* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	matrix_L->multiply_with_vector(u_L, dt_u_L);
	
	for(int i=0; i<ndof;i++){
						lumped= lumped_matrix->get_Ax()[i]*time_step;
						dt_u_L[i] += (vec_bdry->get(i) + vec_source->get(i));
						dt_u_L[i] *=1./lumped;	
 						P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.; flux_scalar[i]=0.;
		}

		//Berechnung von P&Q für gas phase
		for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
					int i = Ai_mass[indx];	
					if(((mass=Ax_mass[indx]*time_step)!=0.)){
						diff=diffusion->get(i,j);
				lumped= lumped_matrix->get_Ax()[i];
	f_rho = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);
				f_rho_v_x = mass*(dt_u_L[i+dof_rho]- dt_u_L[j+dof_rho]) + diff*(u_L[i+dof_rho]- u_L[j+dof_rho]);	
				f_rho_v_y = mass*(dt_u_L[i+dof_v_y_start]- dt_u_L[j+dof_v_y_start]) + diff*(u_L[i+dof_v_y_start]- u_L[j+dof_v_y_start]);	
				f_rhoE = mass*(dt_u_L[i+dof_rhoE_start]- dt_u_L[j+dof_rhoE_start])+ 
							diff*(u_L[i+dof_rhoE_start]- u_L[j+dof_rhoE_start]);			
							
							
		//rho			
				if(f_rho>0.0)	{ 
					P_plus[i]+=f_rho;
				}else if(f_rho<0.0){
				 	P_minus[i]+=f_rho;
				}										
			f = lumped*(u_L[j]-u_L[i]); 
				if(f>Q_plus[i]) Q_plus[i] = f;				
				if(f<Q_minus[i]) Q_minus[i] = f;	
		//v_x
			v_x_i = (u_L[i+dof_rho]/u_L[i]);
				f_v_x= (f_rho_v_x- f_rho*v_x_i)/u_L[i];
				if(f_v_x>0.0)	{ 
					P_plus[i+dof_rho]+=f_v_x;
				}else if(f_v_x<0.0){
				 	P_minus[i+dof_rho]+=f_v_x;
				}										
				f = lumped*((u_L[j+dof_rho]/u_L[j])-v_x_i); 
				if(f>Q_plus[i+dof_rho]) Q_plus[i+dof_rho] = f;				
				if(f<Q_minus[i+dof_rho]) Q_minus[i+dof_rho] = f;	

		//v_y
				v_y_i = (u_L[i+dof_v_y_start]/u_L[i]);
				f_v_y= (f_rho_v_y- f_rho*v_y_i)/u_L[i];
				if(f_v_y>0.0)	{ 
					P_plus[i+dof_v_y_start]+=f_v_y;
				}else if(f_v_x<0.0){
				 	P_minus[i+dof_v_y_start]+=f_v_y;
				}										
				f = lumped*((u_L[j+dof_v_y_start]/u_L[j])-v_y_i); 
				if(f>Q_plus[i+dof_v_y_start]) Q_plus[i+dof_v_y_start] = f;				
				if(f<Q_minus[i+dof_v_y_start]) Q_minus[i+dof_v_y_start] = f;	
		//p
				p_i= (gamma-1)*(u_L[i+dof_rhoE_start]-((v_x_i*v_x_i+v_y_i*v_y_i)*0.5*u_L[i]));
				p_j= (gamma-1)*(u_L[j+dof_rhoE_start]-((u_L[j+dof_rho]*u_L[j+dof_rho]+u_L[j+dof_v_y_start]*u_L[j+dof_v_y_start])/(2*u_L[j])));
				f_p= (gamma-1)*(f_rhoE+ f_rho*(v_x_i*v_x_i+v_y_i*v_y_i)*0.5-(v_x_i*f_rho_v_x+v_y_i*f_rho_v_y));
				if(f_p>0.0)	{ 
					P_plus[i+dof_rhoE_start]+=f_p;
				}else if(f_p<0.0){
				 	P_minus[i+dof_rhoE_start]+=f_p;
				}		
		
				f = lumped*(p_j-p_i); 
				if(f>Q_plus[i+dof_rhoE_start]) Q_plus[i+dof_rhoE_start] = f;				
				if(f<Q_minus[i+dof_rhoE_start]) Q_minus[i+dof_rhoE_start] = f;

					}
				}			
		}

		//Berechnung von P&Q für particle
		for(int j = ndof_g; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx]*time_step)!=0.)&&(j>i)){
								diff = diffusion->get(i,j);				
									f = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);	
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
										if(f>0.0)	{ 
											P_plus[i]+=f;
											P_minus[j]-=f;
										}else if (f<0.0){
										 	P_minus[i]+=f;
											P_plus[j]-=f;
										}
										lumped= lumped_matrix->get_Ax()[i];
										f = lumped*(u_L[j]-u_L[i]); 
										if(f>Q_plus[i]) Q_plus[i] = f;				
										if(f<Q_minus[i]) Q_minus[i] = f;			
										f= lumped*(u_L[i]-u_L[j]);
										if(f>Q_plus[j]) Q_plus[j] = f;	
										if(f<Q_minus[j]) Q_minus[j] = f;

							}	
							
				}			
		}



			//Berechnung von R	
	for(int i=0; i<ndof;i++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
		if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
		if(plus>=1.0) R_plus[i]= 1.0;
		else 	     R_plus[i]= plus;
		if(minus>=1.0) R_minus[i]= 1.0;
		else 	     R_minus[i]= minus;

	}	

	//Berechnung von alpha & f_i
	alpha = 1.0;
		for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx]*time_step)!=0.)&&(i<j)){
								diff = diffusion->get(i,j);
				f_rho_ij = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);
				f_rho_ji = -f_rho_ij;  //antisymmetrisch
				f_rho_v_x = mass*(dt_u_L[i+dof_rho]- dt_u_L[j+dof_rho]) + diff*(u_L[i+dof_rho]- u_L[j+dof_rho]);	//antisymmetrisch
				f_rho_v_y = mass*(dt_u_L[i+dof_v_y_start]- dt_u_L[j+dof_v_y_start]) + diff*(u_L[i+dof_v_y_start]- u_L[j+dof_v_y_start]);//antisymmetrisch	
				f_rhoE = mass*(dt_u_L[i+dof_rhoE_start]- dt_u_L[j+dof_rhoE_start]) + diff*(u_L[i+dof_rhoE_start]- u_L[j+dof_rhoE_start]);	//antisymmetrisch
				
				
				v_x_i = (u_L[i+dof_rho]/u_L[i]);	//vx_i
				v_y_i = (u_L[i+dof_v_y_start]/u_L[i]);
				f_p_ij=(gamma-1)*(f_rhoE+ f_rho*(v_x_i*v_x_i+v_y_i*v_y_i)*0.5-(v_x_i*f_rho_v_x+v_y_i*f_rho_v_y));
				f_v_x_ij = (f_rho_v_x-v_x_i*f_rho_ij)/u_L[i];
				f_v_y_ij = (f_rho_v_y-v_y_i*f_rho_ij)/u_L[i];
				
				v_x_j = (u_L[j+dof_rho]/u_L[j]); //vx_j	
				v_y_j = (u_L[j+dof_v_y_start]/u_L[j]);
				f_p_ji = (gamma-1)*(-f_rhoE +f_rho_ji*(v_x_j*v_x_j+v_y_j*v_y_j)*0.5-(v_x_j*(-f_rho_v_x)+v_y_j*(-f_rho_v_y)));
				f_v_x_ji = (-f_rho_v_x-v_x_j*f_rho_ji)/u_L[j];
				f_v_y_ji = (-f_rho_v_y-v_y_j*f_rho_ji)/u_L[j];
	//density- rho
				if(f_rho_ij>=0.){//R_plus[i]					
						if(f_rho_ji>=0.){//R_plus[j]	
								if(R_plus[i]>R_plus[j]) alpha_rho = R_plus[j];
								else alpha_rho = R_plus[i];
						}else{//R_minus[j]
							if(R_plus[i]>R_minus[j]) alpha_rho = R_minus[j];
							else alpha_rho = R_plus[i];
						}	
				}else{//R_minus[i]
						if(f_rho_ji>=0.){//R_plus[j]	
								if(R_minus[i]>R_plus[j]) alpha_rho = R_plus[j];
								else alpha_rho = R_minus[i];
						}else{//R_minus[j]
							if(R_minus[i]>R_minus[j]) alpha_rho = R_minus[j];
							else alpha_rho = R_minus[i];
						}	
				}		
		//pressure 
				if(f_p_ij>=0.){//R_plus[i]					
						if(f_p_ji>=0.){//R_plus[j]	
								if(R_plus[i+dof_rhoE_start]>R_plus[j+dof_rhoE_start]) alpha_p = R_plus[j+dof_rhoE_start];
								else alpha_p = R_plus[i+dof_rhoE_start];
						}else{//R_minus[j]
							if(R_plus[i+dof_rhoE_start]>R_minus[j+dof_rhoE_start]) alpha_p = R_minus[j+dof_rhoE_start];
							else alpha_p = R_plus[i+dof_rhoE_start];
						}	
				}else{//R_minus[i]
						if(f_p_ji>=0.){//R_plus[j]	
								if(R_minus[i+dof_rhoE_start]>R_plus[j+dof_rhoE_start]) alpha_p = R_plus[j+dof_rhoE_start];
								else alpha_p = R_minus[i+dof_rhoE_start];
						}else{//R_minus[j]
							if(R_minus[i+dof_rhoE_start]>R_minus[j+dof_rhoE_start]) alpha_p = R_minus[j+dof_rhoE_start];
							else alpha_p = R_minus[i+dof_rhoE_start];
						}	
				}
			if(alpha_p>alpha_rho) alpha = alpha_rho;
				else alpha = alpha_p;


				flux_scalar[i] += alpha*f_rho_ij;
				flux_scalar[i+dof_rho ] += alpha*f_rho_v_x;			
				flux_scalar[i+dof_v_y_start ] += alpha*f_rho_v_y;	
				flux_scalar[i+dof_rhoE_start] += alpha*f_rhoE;
				//alpha_ij=alpha_ji & antisymmetric fluxes!		
				flux_scalar[j] -= alpha*f_rho_ij;
				flux_scalar[j+dof_rho ] -= alpha*f_rho_v_x;			
				flux_scalar[j+dof_v_y_start ] -= alpha*f_rho_v_y;	
				flux_scalar[j+dof_rhoE_start] -= alpha*f_rhoE;						
										
							}
				}			
		}


//particle flux
		for(int j = ndof_g; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx]*time_step)!=0.)&&(j>i)){
														diff = diffusion->get(i,j);
								//f = (mass/time_step+ diff/2.)*(u_high[i]- u_high[j])
													//	-(mass/time_step - diff/2.) *(u_old[i]- u_old[j]);
									f = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);	
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
										
										if(f>0.){					
											if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
											else 	alpha = R_plus[i];
										}else if(f<0.){
											if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
											else 	alpha = R_minus[i]; 
										}									
											flux_scalar[i] += alpha*f;
											flux_scalar[j] -= alpha*f;	
									
										
							}
				}			
		}

delete [] dt_u_L ;
}

