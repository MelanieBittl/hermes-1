//Assemble antidiffusive fluxes & Limiter
void antidiffusiveFlux(CSCMatrix<double>* mass_matrix,CSCMatrix<double>* lumped_matrix,CSCMatrix<double>* conv_matrix,CSCMatrix<double>* diffusion,double* u_high, double* u_L, double* u_old,double* flux_scalar, double* P_plus, double* P_minus, double* Q_plus, double* Q_minus,  double* R_plus, double* R_minus, int* smooth_dof=NULL  )
{ 
	int ndof = conv_matrix->get_size();
	double alpha,f, plus, minus,mass, diff;
	double* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;flux_scalar[i]=0.0;}

		//Berechnung von P&Q
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
								diff = diffusion->get(i,j);
								f = (mass/time_step+ diff/2.)*(u_high[i]- u_high[j])
														-(mass/time_step - diff/2.) *(u_old[i]- u_old[j]);	
										if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
										if(f>0.0)	{ 
											P_plus[i]+=f;
											P_minus[j]-=f;
										}else if (f<0.0){
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
														diff = diffusion->get(i,j);
						f = (mass/time_step+ diff/2.)*(u_high[i]- u_high[j])
												-(mass/time_step - diff/2.) *(u_old[i]- u_old[j]);	
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


}

//FCT for lumped projection
template<typename Scalar>
void lumped_flux_limiter(CSCMatrix<Scalar>* mass_matrix,CSCMatrix<Scalar>* lumped_matrix, Scalar* u_L, Scalar* u_H, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus, Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL)
{	int ndof = mass_matrix->get_size();
	Scalar* rhs = new Scalar[ndof];	
	Scalar alpha,f, plus, minus,mass;

	Scalar* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;rhs[i]=0.;}
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

	for(int i=0; i<ndof;i++) u_L[i] += rhs[i]/lumped_matrix->get_Ax()[i];

	delete [] rhs;
}

//MIT FCT_BOOL



//Assemble antidiffusive fluxes & Limiter
void antidiffusiveFlux(bool* fct,CSCMatrix<double>* mass_matrix,CSCMatrix<double>* lumped_matrix,CSCMatrix<double>* conv_matrix,CSCMatrix<double>* diffusion,double* u_high, double* u_L, double* u_old,double* flux_scalar, double* P_plus, double* P_minus, double* Q_plus, double* Q_minus,  double* R_plus, double* R_minus, int* smooth_dof=NULL  )
{ 
	int ndof = conv_matrix->get_size();
	double alpha,f, plus, minus,mass, diff;
	double* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.;flux_scalar[i]=0.0;}

		//Berechnung von P&Q
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				if(fct[j]== false) continue;
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(fct[i]== false) continue;
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
								diff = diffusion->get(i,j);
								f = (mass/time_step+ diff/2.)*(u_high[i]- u_high[j])
														-(mass/time_step - diff/2.) *(u_old[i]- u_old[j]);	
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
		}

			//Berechnung von R	
	for(int i=0; i<ndof;i++){
		if(fct[i]== false) continue;
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
				if(fct[j]== false) continue;
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(fct[i]== false) continue;
							if(((mass=Ax_mass[indx])!=0.)&&(j<i)){
														diff = diffusion->get(i,j);
						f = (mass/time_step+ diff/2.)*(u_high[i]- u_high[j])
												-(mass/time_step - diff/2.) *(u_old[i]- u_old[j]);	
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


}

//FCT for lumped projection
template<typename Scalar>
void lumped_flux_limiter(bool* fct,CSCMatrix<Scalar>* mass_matrix,CSCMatrix<Scalar>* lumped_matrix, Scalar* u_L, Scalar* u_H, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,  Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL)
{	
	int ndof = mass_matrix->get_size();
	Scalar* rhs = new Scalar[ndof];
	lumped_matrix->multiply_with_vector(u_L, rhs); 

	Scalar alpha,f, plus, minus,mass;

	Scalar* Ax_mass = mass_matrix->get_Ax();
	int* Ai_mass = mass_matrix->get_Ai();
	int* Ap_mass = mass_matrix->get_Ap();

	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.;}
		//Berechnung von P&Q
		for(int j = 0; j<ndof; j++){ //Spalten durchlaufen
				if(fct[j]== false) continue;
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];
							if(fct[i]== false) continue;	
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
										f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i]); 
										if(f>Q_plus[i]) Q_plus[i] = f;				
										if(f<Q_minus[i]) Q_minus[i] = f;			
										f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j]);
										if(f>Q_plus[j]) Q_plus[j] = f;	
										if(f<Q_minus[j]) Q_minus[j] = f;	
							}
				}			
		}


		//Berechnung von R
		for(int i=0; i<ndof;i++){
			if(fct[i]== false) continue;	
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
				if(fct[j]== false) continue;
				for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++){	
							int i = Ai_mass[indx];	
							if(fct[i]== false) continue;
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



	SimpleVector<Scalar>* vec_rhs = new SimpleVector<Scalar>(ndof);	
	vec_rhs->zero(); vec_rhs->add_vector(rhs);
	Scalar* sol =NULL;
	UMFPackLinearMatrixSolver<Scalar>* lowOrd = new UMFPackLinearMatrixSolver<Scalar>(lumped_matrix,vec_rhs);
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
	delete vec_rhs;
	delete [] rhs;
}










