
//Assemble antidiffusive fluxes & Limiter
template<typename Scalar>
void antidiffusiveFlux(UMFPackMatrix<Scalar>* mass_matrix,UMFPackMatrix<Scalar>* lumped_matrix,UMFPackMatrix<Scalar>* conv_matrix,UMFPackMatrix<Scalar>* diffusion,UMFPackVector<Scalar>* flux_dt_rhs, Scalar* u_L, Scalar* u_old,
Scalar* flux_scalar, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,Scalar* Q_plus_old, Scalar* Q_minus_old, Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL  )
{ 


	int ndof = conv_matrix->get_size();
	Scalar flux_dt_scalar[ndof];	
	LinearSolver<Scalar>* flux_dt;
	Scalar* dt_u_L = NULL;
	Scalar alpha,f, plus, minus,mass;

	Scalar* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++) flux_scalar[i]=0.0;

			conv_matrix->multiply_with_vector(u_L, flux_dt_scalar);
			flux_dt_rhs->zero(); flux_dt_rhs->add_vector(flux_dt_scalar);  //K u^L	
			flux_dt = create_linear_solver(matrix_solver,mass_matrix,flux_dt_rhs); //M_c u_t = K u^L
			if(flux_dt->solve())	dt_u_L = flux_dt->get_sln_vector();	
				else error ("Matrix solver failed.\n");
	
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=Q_plus_old[i];Q_minus[i]=Q_minus_old[i];}

		//Berechnung von P&Q
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
										f = mass*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
										
											if(f>0.0)	{ 
												P_plus[i]+=f;
												P_minus[j]-=f;
											}else if(f<0.0){
											 	P_minus[i]+=f;
												P_plus[j]-=f;
											}
										
										f = lumped_matrix->get_Ax()[i]*(u_L[j]-u_L[i])/time_step; 
										if(f>Q_plus[i]) Q_plus[i] = f;				
										if(f<Q_minus[i]) Q_minus[i] = f;			
										f= lumped_matrix->get_Ax()[j]*(u_L[i]-u_L[j])/time_step;
										if(f>Q_plus[j]) Q_plus[j] = f;	
										if(f<Q_minus[j]) Q_minus[j] = f;	
							}
				}			
		}
		//Berechnung von P&Q
		/*for(int i =0;i<ndof;i++){
			for(int j =(i+1); j<ndof;j++){		
				if((mass=mass_matrix->get(i,j))!=0.0){		
						f = mass*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
					if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0)	{ 
						P_plus[i]+=f;
						P_minus[j]-=f;
					}else if (f<0.0){
					 	P_minus[i]+=f;
						P_plus[j]-=f;
					}
					f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i])/time_step; 
					if(f>Q_plus[i]) Q_plus[i] = f;				
					if(f<Q_minus[i]) Q_minus[i] = f;			
					f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j])/time_step;
					if(f>Q_plus[j]) Q_plus[j] = f;	
					if(f<Q_minus[j]) Q_minus[j] = f;	

				}
			}
		}*/

			//Berechnung von R	
	for(int i=0; i<ndof;i++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
		if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
		if(plus>=1.0) R_plus[i]= 1.0;
		else 	     R_plus[i]= plus;
		if(minus>=1.0) R_minus[i]= 1.0;
		else 	     R_minus[i]= minus;
		if(smooth_dof!=NULL){
				if(smooth_dof[i]==1){ R_plus[i]= 1.0;R_minus[i]= 1.0;}
		}	
	}	

	//Berechnung von alpha & f_i
	alpha = 1.0;
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
										f = mass*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
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


/*
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){	
					if((mass=mass_matrix->get(i,j))!=0.0){		
						f = mass*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
						if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
						if(f>0){					
							if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
							else 	alpha = R_plus[i];
						}else{
							if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
							else 	alpha = R_minus[i]; 
						}
						flux_scalar[i] += alpha*f;
						flux_scalar[j] -= alpha*f;
					}				
				}
			}
*/



			//Cleanup	
			delete flux_dt;
		

}


//FCT for lumped projection
template<typename Scalar>
void lumped_flux_limiter(UMFPackMatrix<Scalar>* mass_matrix,UMFPackMatrix<Scalar>* lumped_matrix, Scalar* u_L, Scalar* u_H, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,Scalar* Q_plus_old, Scalar* Q_minus_old,  Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL)
{	int ndof = mass_matrix->get_size();
	UMFPackVector<Scalar>* vec_rhs = new UMFPackVector<Scalar>(ndof);	
	vec_rhs->zero(); 
	Scalar* rhs = new Scalar[ndof];
	lumped_matrix->multiply_with_vector(u_L, rhs); 
	Scalar alpha,f, plus, minus,mass;

	Scalar* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=Q_plus_old[i]*time_step;Q_minus[i]=Q_minus_old[i]*time_step;}

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


/*
	for(int i =0;i<ndof;i++){
			for(int j =(i+1); j<ndof;j++){	
				if((mass=mass_matrix->get(i,j))!=0.0){	
						f = mass*(u_H[i]- u_H[j]);
						if( (f*(u_H[j]- u_H[i])) > 0.0) f = 0.0; //prelimiting step
						if(f>0.0)	{ 
							P_plus[i]+=f;
							P_minus[j]-=f;
						}else if (f<0.0){
						 	P_minus[i]+=f;
							P_plus[j]-=f;
						}
						f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i]); 
						if(f>Q_plus[i]) Q_plus[i] = f;				
						if(f<Q_minus[i]) Q_minus[i] = f;			
						f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j]); 
						if(f>Q_plus[j]) Q_plus[j] = f;	
						if(f<Q_minus[j]) Q_minus[j] = f;		
				}	
			}
		}
*/

		//Berechnung von R
		for(int i=0; i<ndof;i++){
			plus = 1.0; minus = 1.0;		
			if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
			if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
			if(plus>=1.0) R_plus[i]= 1.0;
			else 	     R_plus[i]= plus;
			if(minus>=1.0) R_minus[i]= 1.0;
			else 	     R_minus[i]= minus;
		if(smooth_dof!=NULL){
				if(smooth_dof[i]==1){ R_plus[i]= 1.0;R_minus[i]= 1.0;}
		}			
		}	

	//Berechnung von alpha & f_i
	alpha = 1.0;
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
										f = mass*(u_H[i]- u_H[j]);
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




/*	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){	
				if((mass=mass_matrix->get(i,j))!=0.0){			
					f = mass*(u_H[i]- u_H[j]);		
					if( (f*(u_H[j]- u_H[i])) > 0.0) f = 0.0; //prelimiting step
					if(f>=0.0){					
						if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
						else 	alpha = R_plus[i];
					}else {
						if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
						else 	alpha = R_minus[i]; 
					}
					rhs[i]+= alpha*f;
					rhs[j]-= alpha*f;
				}			
		}
	}*/




	vec_rhs->add_vector(rhs);
	double* sol =NULL;
	UMFPackLinearSolver<double>* lowOrd = new UMFPackLinearSolver<double>(lumped_matrix,vec_rhs);
	if(lowOrd->solve()){ 
		sol = lowOrd->get_sln_vector();  			
	}else error ("Matrix in lumped_flux solver failed.\n");
	for(int i=0; i<ndof;i++) u_L[i] =sol[i];

	delete lowOrd;
	
	delete vec_rhs;
	delete [] rhs;
}

