
//Assemble antidiffusive fluxes & Limiter
//f_ij und alpha_ij werden nicht explizit berechnet!! da scalar** flux = new_matrix<scalar>(ndof,ndof); zuviel Speicher braucht
template<typename Scalar>
void antidiffusiveFlux(AsmList<Scalar>* al,AsmList<Scalar>* al_2,UMFPackMatrix<Scalar>* mass_matrix,UMFPackMatrix<Scalar>* lumped_matrix,UMFPackMatrix<Scalar>* conv_matrix,UMFPackMatrix<Scalar>* diffusion,UMFPackVector<Scalar>* flux_dt_rhs, Scalar* u_L, Scalar* flux_scalar, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,Scalar* Q_plus_old, Scalar* Q_minus_old, Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL  )
{ //al==NULL =>flux=0
	int ndof = conv_matrix->get_size();
	Scalar flux_dt_scalar[ndof];	
	LinearSolver<Scalar>* flux_dt;
	Scalar* dt_u_L = NULL;
	Scalar alpha,f, plus, minus,mass;
	for(int i=0; i<ndof;i++) flux_scalar[i]=0.0;
	if(al!=NULL){
			conv_matrix->multiply_with_vector(u_L, flux_dt_scalar);
			flux_dt_rhs->zero(); flux_dt_rhs->add_vector(flux_dt_scalar);  //K u^L	
			flux_dt = create_linear_solver(matrix_solver,mass_matrix,flux_dt_rhs); //M_c u_t = K u^L
			if(flux_dt->solve())	dt_u_L = flux_dt->get_sln_vector();	
				else error ("Matrix solver failed.\n");

			//Berechnung von P&Q
			for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=Q_plus_old[i];Q_minus[i]=Q_minus_old[i];}
			for(unsigned int i = 0; i < al->get_cnt(); i ++){
			 	for(unsigned int j = (i+1); j < al->get_cnt(); j ++){
					if((al->get_dof()[i]!=al->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al->get_dof()[j]))!=0.0)){		
						f = mass*(dt_u_L[al->get_dof()[i]]- dt_u_L[al->get_dof()[j]]) 
							+ diffusion->get(al->get_dof()[i],al->get_dof()[j])*(u_L[al->get_dof()[i]]- u_L[al->get_dof()[j]]);
						if( (f*(u_L[al->get_dof()[j]]- u_L[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
						if(f>0.0)	{ 
							P_plus[al->get_dof()[i]]+=f;
							P_minus[al->get_dof()[j]]-=f;
						}else if (f<0.0){
						 	P_minus[al->get_dof()[i]]+=f;
							P_plus[al->get_dof()[j]]-=f;
						}
						f = lumped_matrix->get(al->get_dof()[i],al->get_dof()[i])*(u_L[al->get_dof()[j]]-u_L[al->get_dof()[i]])/time_step; 
						if(f>Q_plus[al->get_dof()[i]]) Q_plus[al->get_dof()[i]] = f;				
						if(f<Q_minus[al->get_dof()[i]]) Q_minus[al->get_dof()[i]] = f;			
						f= lumped_matrix->get(al->get_dof()[j],al->get_dof()[j])*(u_L[al->get_dof()[i]]-u_L[al->get_dof()[j]])/time_step; 
						if(f>Q_plus[al->get_dof()[j]]) Q_plus[al->get_dof()[j]] = f;	
						if(f<Q_minus[al->get_dof()[j]]) Q_minus[al->get_dof()[j]] = f;
					}
				}
			}
		//fuer list_2
			for(unsigned int i = 0; i < al->get_cnt(); i ++){
			 	for(unsigned int j = 0; j < al_2->get_cnt(); j ++){
					if((al->get_dof()[i]!=al_2->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al_2->get_dof()[j]))!=0.0)){			
						f = mass*(dt_u_L[al->get_dof()[i]]- dt_u_L[al_2->get_dof()[j]]) 
							+ diffusion->get(al->get_dof()[i],al_2->get_dof()[j])*(u_L[al->get_dof()[i]]- u_L[al_2->get_dof()[j]]);
						if( (f*(u_L[al_2->get_dof()[j]]- u_L[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
						if(f>0.0)	{ 
							P_plus[al->get_dof()[i]]+=f;
							P_minus[al_2->get_dof()[j]]-=f;
						}else if (f<0.0){
						 	P_minus[al->get_dof()[i]]+=f;
							P_plus[al_2->get_dof()[j]]-=f;
						}
						f = lumped_matrix->get(al->get_dof()[i],al->get_dof()[i])*(u_L[al_2->get_dof()[j]]-u_L[al->get_dof()[i]])/time_step; 
						if(f>Q_plus[al->get_dof()[i]]) Q_plus[al->get_dof()[i]] = f;				
						if(f<Q_minus[al->get_dof()[i]]) Q_minus[al->get_dof()[i]] = f;			
						f= lumped_matrix->get(al_2->get_dof()[j],al_2->get_dof()[j])*(u_L[al->get_dof()[i]]-u_L[al_2->get_dof()[j]])/time_step; 
						if(f>Q_plus[al_2->get_dof()[j]]) Q_plus[al_2->get_dof()[j]] = f;	
						if(f<Q_minus[al_2->get_dof()[j]]) Q_minus[al_2->get_dof()[j]] = f;
					}
				}
			}
			



			//Berechnung von R	
			for(unsigned int i = 0; i < al->get_cnt(); i ++){
				plus = 1.0; minus = 1.0;		
				if(P_plus[al->get_dof()[i]]!=0.0)  plus = Q_plus[al->get_dof()[i]]/P_plus[al->get_dof()[i]];		
				if(P_minus[al->get_dof()[i]]!=0.0) minus = Q_minus[al->get_dof()[i]]/P_minus[al->get_dof()[i]];			
				if(plus>=1.0) R_plus[al->get_dof()[i]]= 1.0;
				else 	     R_plus[al->get_dof()[i]]= plus;
				if(minus>=1.0) R_minus[al->get_dof()[i]]= 1.0;
				else 	     R_minus[al->get_dof()[i]]= minus;

				if(smooth_dof!=NULL){
						if(smooth_dof[al->get_dof()[i]]==1){ R_plus[al->get_dof()[i]]= 1.0;R_minus[al->get_dof()[i]]= 1.0;}
				}	
			}

			//fuer list2
			for(unsigned int i = 0; i < al_2->get_cnt(); i ++){
				plus = 1.0; minus = 1.0;		
				/*if(P_plus[al_2->get_dof()[i]]!=0.0)  plus = Q_plus[al_2->get_dof()[i]]/P_plus[al_2->get_dof()[i]];		
				if(P_minus[al_2->get_dof()[i]]!=0.0) minus = Q_minus[al_2->get_dof()[i]]/P_minus[al_2->get_dof()[i]];			
				if(plus>=1.0) R_plus[al_2->get_dof()[i]]= 1.0;
				else 	     R_plus[al_2->get_dof()[i]]= plus;
				if(minus>=1.0) R_minus[al_2->get_dof()[i]]= 1.0;
				else 	     R_minus[al_2->get_dof()[i]]= minus;	*/
					R_plus[al_2->get_dof()[i]]= 1.0;
				R_minus[al_2->get_dof()[i]]= 1.0;
				//if(smooth_dof!=NULL){
				//		if(smooth_dof[al_2->get_dof()[i]]==1){ R_plus[al_2->get_dof()[i]]= 1.0;R_minus[al_2->get_dof()[i]]= 1.0;}
				//}	
	
			}

	
			//Berechnung von alpha & f_i
			alpha = 1.0;
			for(unsigned int i = 0; i < al->get_cnt(); i ++){
			 	for(unsigned int j = (i+1); j < al->get_cnt(); j ++){	
					if((al->get_dof()[i]!=al->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al->get_dof()[j]))!=0.0)){		
					f= mass*(dt_u_L[al->get_dof()[i]]- dt_u_L[al->get_dof()[j]]) 
						+ diffusion->get(al->get_dof()[i],al->get_dof()[j])*(u_L[al->get_dof()[i]]- u_L[al->get_dof()[j]]);	
					if( (f*(u_L[al->get_dof()[j]]- u_L[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0){					
						if(R_plus[al->get_dof()[i]]>R_minus[al->get_dof()[j]]) alpha = R_minus[al->get_dof()[j]];
						else 	alpha = R_plus[al->get_dof()[i]];
					}else if (f<0.0){
						if(R_minus[al->get_dof()[i]]>R_plus[al->get_dof()[j]]) alpha = R_plus[al->get_dof()[j]];
						else 	alpha = R_minus[al->get_dof()[i]]; 
					}
					flux_scalar[al->get_dof()[i]] += alpha*f;
					flux_scalar[al->get_dof()[j]] -= alpha*f;
					}				
				}
			}
			//fuer list2
			for(unsigned int i = 0; i < al->get_cnt(); i ++){
			 	for(unsigned int j = 0; j < al_2->get_cnt(); j ++){	
					if((al->get_dof()[i]!=al_2->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al_2->get_dof()[j]))!=0.0)){				
					f= mass*(dt_u_L[al->get_dof()[i]]- dt_u_L[al_2->get_dof()[j]]) 
						+ diffusion->get(al->get_dof()[i],al_2->get_dof()[j])*(u_L[al->get_dof()[i]]- u_L[al_2->get_dof()[j]]);	
					if( (f*(u_L[al_2->get_dof()[j]]- u_L[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0){					
						if(R_plus[al->get_dof()[i]]>R_minus[al_2->get_dof()[j]]) alpha = R_minus[al_2->get_dof()[j]];
						else 	alpha = R_plus[al->get_dof()[i]];
					}else if (f<0.0){
						if(R_minus[al->get_dof()[i]]>R_plus[al_2->get_dof()[j]]) alpha = R_plus[al_2->get_dof()[j]];
						else 	alpha = R_minus[al->get_dof()[i]]; 
					}
					flux_scalar[al->get_dof()[i]] += alpha*f;     
					flux_scalar[al_2->get_dof()[j]] -= alpha*f;  
					}				
				}
			}



			//Cleanup	
			delete flux_dt;
	}	

}


//FCT for lumped projection
template<typename Scalar>
void lumped_flux_limiter(AsmList<Scalar>* al,AsmList<Scalar>* al_2,UMFPackMatrix<Scalar>* mass_matrix,UMFPackMatrix<Scalar>* lumped_matrix, Scalar* u_L, Scalar* u_H, Scalar* P_plus, Scalar* P_minus, Scalar* Q_plus, Scalar* Q_minus,Scalar* Q_plus_old, Scalar* Q_minus_old,  Scalar* R_plus, Scalar* R_minus, int* smooth_dof=NULL)
{
	int ndof = mass_matrix->get_size();
	UMFPackVector<Scalar>* vec_rhs = new UMFPackVector<Scalar>(ndof);	
	vec_rhs->zero(); 
	Scalar* rhs = new Scalar[ndof];
	lumped_matrix->multiply_with_vector(u_L, rhs); 

	Scalar alpha,f, plus, minus,mass;
	if(al!=NULL){
		//Berechnung von P&Q
			for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=Q_plus_old[i]*time_step;Q_minus[i]=Q_minus_old[i]*time_step;}
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
		 	for(unsigned int j = (i+1); j < al->get_cnt(); j ++){
					if((al->get_dof()[i]!=al->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al->get_dof()[j]))!=0.0)){	
					f = mass*(u_H[al->get_dof()[i]]- u_H[al->get_dof()[j]]);								
					if( (f*(u_H[al->get_dof()[j]]- u_H[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0)	{ 
						P_plus[al->get_dof()[i]]+=f;
						P_minus[al->get_dof()[j]]-=f;
					}else if (f<0.0){
					 	P_minus[al->get_dof()[i]]+=f;
						P_plus[al->get_dof()[j]]-=f;
					}				
					f = lumped_matrix->get(al->get_dof()[i],al->get_dof()[i])*(u_L[al->get_dof()[j]]-u_L[al->get_dof()[i]]); 
					if(f>Q_plus[al->get_dof()[i]]) Q_plus[al->get_dof()[i]] = f;				
					if(f<Q_minus[al->get_dof()[i]]) Q_minus[al->get_dof()[i]] = f;
					f= lumped_matrix->get(al->get_dof()[j],al->get_dof()[j])*(u_L[al->get_dof()[i]]-u_L[al->get_dof()[j]]);
					if(f>Q_plus[al->get_dof()[j]]) Q_plus[al->get_dof()[j]] = f;	
					if(f<Q_minus[al->get_dof()[j]]) Q_minus[al->get_dof()[j]] = f;
				}
			}
		}
	//fuer list2
	for(unsigned int i = 0; i < al->get_cnt(); i ++){
		 	for(unsigned int j = 0; j < al_2->get_cnt(); j ++){
			if((al->get_dof()[i]!=al_2->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al_2->get_dof()[j]))!=0.0)){		
					f = mass*(u_H[al->get_dof()[i]]- u_H[al_2->get_dof()[j]]);				
					if( (f*(u_H[al_2->get_dof()[j]]- u_H[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0)	{ 
						P_plus[al->get_dof()[i]]+=f;
						P_minus[al_2->get_dof()[j]]-=f;
					}else if (f<0.0){
					 	P_minus[al->get_dof()[i]]+=f;
						P_plus[al_2->get_dof()[j]]-=f;
					}				
					f = lumped_matrix->get(al->get_dof()[i],al->get_dof()[i])*(u_L[al_2->get_dof()[j]]-u_L[al->get_dof()[i]]); 
					if(f>Q_plus[al->get_dof()[i]]) Q_plus[al->get_dof()[i]] = f;				
					if(f<Q_minus[al->get_dof()[i]]) Q_minus[al->get_dof()[i]] = f;
					f= lumped_matrix->get(al->get_dof()[j],al->get_dof()[j])*(u_L[al->get_dof()[i]]-u_L[al_2->get_dof()[j]]);
					if(f>Q_plus[al_2->get_dof()[j]]) Q_plus[al_2->get_dof()[j]] = f;	
					if(f<Q_minus[al_2->get_dof()[j]]) Q_minus[al_2->get_dof()[j]] = f;
				}
			}
		}


		//Berechnung von R
	for(unsigned int i = 0; i < al->get_cnt(); i ++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[al->get_dof()[i]]!=0.0)  plus = Q_plus[al->get_dof()[i]]/P_plus[al->get_dof()[i]];		
		if(P_minus[al->get_dof()[i]]!=0.0) minus =Q_minus[al->get_dof()[i]]/P_minus[al->get_dof()[i]];			
		if(plus>=1.0) R_plus[al->get_dof()[i]]= 1.0;
		else 	     R_plus[al->get_dof()[i]]= plus;
		if(minus>=1.0) R_minus[al->get_dof()[i]]= 1.0;
		else 	     R_minus[al->get_dof()[i]]= minus;	
					if(smooth_dof!=NULL){
						if(smooth_dof[al->get_dof()[i]]==1){ R_plus[al->get_dof()[i]]= 1.0;R_minus[al->get_dof()[i]]= 1.0;}
				}	
	}	
	//fuer list2, R =1 setzen;
	for(unsigned int i = 0; i < al_2->get_cnt(); i ++){
		plus = 1.0; minus = 1.0;		
/*
		if(P_plus[al_2->get_dof()[i]]!=0.0)  plus = lumped_matrix->get(al_2->get_dof()[i],al_2->get_dof()[i])*Q_plus[al_2->get_dof()[i]]/P_plus[al_2->get_dof()[i]];		
		if(P_minus[al_2->get_dof()[i]]!=0.0) minus = lumped_matrix->get(al_2->get_dof()[i],al_2->get_dof()[i])*Q_minus[al_2->get_dof()[i]]/P_minus[al_2->get_dof()[i]];		
		if(plus>=1.0) R_plus[al_2->get_dof()[i]]= 1.0;
		else 	     R_plus[al_2->get_dof()[i]]= plus;
		if(minus>=1.0) R_minus[al_2->get_dof()[i]]= 1.0;
		else 	     R_minus[al_2->get_dof()[i]]= minus;	*/

			R_plus[al_2->get_dof()[i]]= 1.0;
			R_minus[al_2->get_dof()[i]]= 1.0;
	}

		//Berechnung von alpha & f_i
		alpha = 1.0;
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
		 	for(unsigned int j = (i+1); j < al->get_cnt(); j ++){	
					if((al->get_dof()[i]!=al->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al->get_dof()[j]))!=0.0)){			
				f= mass*(u_H[al->get_dof()[i]]- u_H[al->get_dof()[j]]) ;	
				if( (f*(u_H[al->get_dof()[j]]- u_H[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step				
				if(f>0.0){					
					if(R_plus[al->get_dof()[i]]>R_minus[al->get_dof()[j]]) alpha = R_minus[al->get_dof()[j]];
					else 	alpha = R_plus[al->get_dof()[i]];
				}else if(f<0.0){
					if(R_minus[al->get_dof()[i]]>R_plus[al->get_dof()[j]]) alpha = R_plus[al->get_dof()[j]];
					else 	alpha = R_minus[al->get_dof()[i]]; 
				}
					rhs[al->get_dof()[i]]+= alpha*f;
					rhs[al->get_dof()[j]]-= alpha*f;
			  }				
			}
		}


		//fuer list2
				alpha = 1.0;
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
		 	for(unsigned int j = 0; j < al_2->get_cnt(); j ++){	
					if((al->get_dof()[i]!=al_2->get_dof()[j])&&((mass=mass_matrix->get(al->get_dof()[i],al_2->get_dof()[j]))!=0.0)){						
				f= mass*(u_H[al->get_dof()[i]]- u_H[al_2->get_dof()[j]]) ;	
				if( (f*(u_H[al_2->get_dof()[j]]- u_H[al->get_dof()[i]])) > 0.0) f = 0.0; //prelimiting step				
				if(f>0.0){					
					if(R_plus[al->get_dof()[i]]>R_minus[al_2->get_dof()[j]]) alpha = R_minus[al_2->get_dof()[j]];
					else 	alpha = R_plus[al->get_dof()[i]];
				}else if(f<0.0){
					if(R_minus[al->get_dof()[i]]>R_plus[al_2->get_dof()[j]]) alpha = R_plus[al_2->get_dof()[j]];
					else 	alpha = R_minus[al->get_dof()[i]]; 
				}				
					rhs[al->get_dof()[i]]+= alpha*f;
					rhs[al_2->get_dof()[j]]-= alpha*f;


			  }				
			}
		}
		

	vec_rhs->add_vector(rhs);
	Scalar* sol =NULL;
	UMFPackLinearSolver<Scalar>* lowOrd = new UMFPackLinearSolver<Scalar>(lumped_matrix,vec_rhs);
	if(lowOrd->solve()){ 
		sol = lowOrd->get_sln_vector();  			
	}else error ("Matrix in lumped_flux solver failed.\n");
	for(int i=0; i<ndof;i++) u_L[i] =sol[i];

	delete lowOrd;
		
	}

	
	delete vec_rhs;
	delete [] rhs;

}

