
//Assemble antidiffusive fluxes & Limiter
void antidiffusiveFlux(CSCMatrix<double>* mass_matrix,CSCMatrix<double>* lumped_matrix,CSCMatrix<double>* diffusion,CSCMatrix<double>* matrix_L,	double* u_L, double* flux_scalar, double* P_plus, double* P_minus, double* Q_plus, double* Q_minus, double* R_plus, double* R_minus, int dof_rho, int dof_v_x, int dof_v_y, int dof_e )
{ 

	int ndof = diffusion->get_size();
	double* dt_u_L = new double[ndof];
	double alpha,f, plus, minus,mass, diff,lumped;
int dof_rhoE_start = dof_rho + dof_v_x + dof_v_y;
int dof_v_y_start = dof_rho + dof_v_x ;

	double f_rho_ij, f_rho_ji, alpha_rho;
	double f_p_ij, f_p_ji, alpha_p;

	double f_rho, f_rho_x, f_rho_y, f_rhoE;
	double f_x,f_y,f_p, v_x, v_y,p_i,p_j;

	double* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	matrix_L->multiply_with_vector(u_L, dt_u_L);
	for(int i=0; i<ndof;i++){
						dt_u_L[i] *=1./lumped_matrix->get(i,i);	
 						P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.; flux_scalar[i]=0.;
		}

		//Berechnung von P&Q
		for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
					int i = Ai_mass[indx];	
					if(((mass=Ax_mass[indx])!=0.)){
							diff=diffusion->get(i,j);
							lumped= lumped_matrix->get_Ax()[i];
							f_rho = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);
							f_rho_x = mass*(dt_u_L[i+dof_rho]- dt_u_L[j+dof_rho]) + diff*(u_L[i+dof_rho]- u_L[j+dof_rho]);	
							f_rho_y = mass*(dt_u_L[i+dof_v_y_start]- dt_u_L[j+dof_v_y_start]) + diff*(u_L[i+dof_v_y_start]- u_L[j+dof_v_y_start]);	
							f_rhoE = mass*(dt_u_L[i+dof_rhoE_start]- dt_u_L[j+dof_rhoE_start]) + diff*(u_L[i+dof_rhoE_start]- u_L[j+dof_rhoE_start]);			
			//rho			
							if(f_rho>0.0)	{ 
								P_plus[i]+=f_rho;
							}else if(f_rho<0.0){
							 	P_minus[i]+=f_rho;
							}										
						f = lumped*(u_L[j]-u_L[i])/time_step; 
						if(f>Q_plus[i]) Q_plus[i] = f;				
						if(f<Q_minus[i]) Q_minus[i] = f;	
		//v_x
					v_x = (u_L[i+dof_rho]/u_L[i]);
					f_x= (f_rho_x- f_rho*v_x)/u_L[i];
							if(f_x>0.0)	{ 
								P_plus[i+dof_rho]+=f_x;
							}else if(f_x<0.0){
							 	P_minus[i+dof_rho]+=f_x;
							}										
						f = lumped*((u_L[j+dof_rho]/u_L[j])-v_x)/time_step; 
						if(f>Q_plus[i+dof_rho]) Q_plus[i+dof_rho] = f;				
						if(f<Q_minus[i+dof_rho]) Q_minus[i+dof_rho] = f;	

			//v_y
					v_y = (u_L[i+dof_v_y_start]/u_L[i]);
					f_y= (f_rho_y- f_rho*v_y)/u_L[i];
							if(f_y>0.0)	{ 
								P_plus[i+dof_v_y_start]+=f_y;
							}else if(f_x<0.0){
							 	P_minus[i+dof_v_y_start]+=f_y;
							}										
						f = lumped*((u_L[j+dof_v_y_start]/u_L[j])-v_y)/time_step; 
						if(f>Q_plus[i+dof_v_y_start]) Q_plus[i+dof_v_y_start] = f;				
						if(f<Q_minus[i+dof_v_y_start]) Q_minus[i+dof_v_y_start] = f;	
			//p
					p_i= (KAPPA-1)*(u_L[i+dof_rhoE_start]-((v_x*v_x+v_y*v_y)*0.5*u_L[i]));
					p_j= (KAPPA-1)*(u_L[j+dof_rhoE_start]-((u_L[j+dof_rho]*u_L[j+dof_rho]+u_L[j+dof_v_y_start]*u_L[j+dof_v_y_start])/(2*u_L[j])));
					f_p= (KAPPA-1)*(f_rhoE+ f_rho*(v_x*v_x+v_y*v_y)*0.5-(v_x*f_rho_x+v_y*f_rho_y));
							if(f_p>0.0)	{ 
								P_plus[i+dof_rhoE_start]+=f_p;
							}else if(f_p<0.0){
							 	P_minus[i+dof_rhoE_start]+=f_p;
							}										
						f = lumped*(p_j-p_i)/time_step; 
						if(f>Q_plus[i+dof_rhoE_start]) Q_plus[i+dof_rhoE_start] = f;				
						if(f<Q_minus[i+dof_rhoE_start]) Q_minus[i+dof_rhoE_start] = f;	


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
							if(((mass=Ax_mass[indx])!=0.)&&(j<dof_rho)){
										diff = diffusion->get(i,j);
							f_rho_ij = mass*(dt_u_L[i]- dt_u_L[j]) + diff*(u_L[i]- u_L[j]);
							f_rho_ji = -f_rho_ij;  //antisymmetrisch
							f_rho_x = mass*(dt_u_L[i+dof_rho]- dt_u_L[j+dof_rho]) + diff*(u_L[i+dof_rho]- u_L[j+dof_rho]);	//antisymmetrisch
							f_rho_y = mass*(dt_u_L[i+dof_v_y_start]- dt_u_L[j+dof_v_y_start]) + diff*(u_L[i+dof_v_y_start]- u_L[j+dof_v_y_start]);//antisymmetrisch	
							f_rhoE = mass*(dt_u_L[i+dof_rhoE_start]- dt_u_L[j+dof_rhoE_start]) + diff*(u_L[i+dof_rhoE_start]- u_L[j+dof_rhoE_start]);	//antisymmetrisch
							v_x = (u_L[i+dof_rho]/u_L[i]);	
							v_y = (u_L[i+dof_v_y_start]/u_L[i]);
							f_p_ij=(KAPPA-1)*(f_rhoE+ f_rho*(v_x*v_x+v_y*v_y)*0.5-(v_x*f_rho_x+v_y*f_rho_y));
							v_x = (u_L[j+dof_rho]/u_L[j]);	
							v_y = (u_L[j+dof_v_y_start]/u_L[j]);
							f_p_ji = (KAPPA-1)*(-f_rhoE +f_rho_ji*(v_x*v_x+v_y*v_y)*0.5-(v_x*(-f_rho_x)+v_y*(-f_rho_y)));
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
										flux_scalar[i+dof_rho ] += alpha*f_rho_x;			
										flux_scalar[i+dof_v_y_start ] += alpha*f_rho_y;	
										flux_scalar[i+dof_rhoE_start] += alpha*f_rhoE;						
										
							}
				}			
		}
delete [] dt_u_L ;

}

//FCT for lumped projection
template<typename Scalar>
void lumped_flux_limiter(CSCMatrix<Scalar>* mass_matrix,CSCMatrix<Scalar>* lumped_matrix, Scalar* u_L, Scalar* u_H, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,Scalar* R_plus, Scalar* R_minus, int dof_rho, int dof_v_x, int dof_v_y, int dof_e )
{	int ndof = mass_matrix->get_size();
	SimpleVector<Scalar>* vec_rhs = new SimpleVector<Scalar>(ndof);	
	vec_rhs->zero(); 
	Scalar* rhs = new Scalar[ndof];
	lumped_matrix->multiply_with_vector(u_L, rhs); 
	Scalar alpha,f, plus, minus,mass;

	Scalar* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.;}


		//Berechnung von P&Q
			for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
										f = mass*(u_H[i]- u_H[j]);
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
										if(f>0.0)	{ 
											P_plus[i]+=f;
											P_minus[j]-=f;
										}else if(f<0.0){
										 	P_minus[i]+=f;
											P_plus[j]-=f;
										}										
										f = lumped_matrix->get_Ax()[i]*(u_L[j]-u_L[i]); 
										if(f>Q_plus[i]) Q_plus[i] = f;				
										if(f<Q_minus[i]) Q_minus[i] = f;			
										f= lumped_matrix->get_Ax()[j]*(u_L[i]-u_L[j]);
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
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
										f = mass*(u_H[i]- u_H[j]);
										if((j<dof_rho)||(j>=(dof_rho+dof_v_x+dof_v_y))){			
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
										if(f>0.){					
											if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
											else 	alpha = R_plus[i];
										}else if(f<0.){
											if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
											else 	alpha = R_minus[i]; 
										}					
												rhs[i] += alpha*f;
												rhs[j] -= alpha*f;
										}									
										
							}
				}			
		}


	vec_rhs->add_vector(rhs);
	double* sol =NULL;
	UMFPackLinearMatrixSolver<double>* lowOrd = new UMFPackLinearMatrixSolver<double>(lumped_matrix,vec_rhs);
  try
  {
   lowOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	 
		sol = lowOrd->get_sln_vector();  			
	
	for(int i=0; i<ndof;i++) u_L[i] =sol[i];
delete lowOrd;



//for(int i=0; i<ndof;i++) u_L[i] = vec_rhs->get(i)/lumped_matrix->get(i,i);

	
	
	delete vec_rhs;
	delete [] rhs;
}

