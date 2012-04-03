
//Mass lumping
UMFPackMatrix<double>* massLumping_simple(UMFPackMatrix<double>* mass_matrix){
	 UMFPackMatrix<double>* lumped_matrix = new UMFPackMatrix<double>;   //M_L
	 int size = mass_matrix->get_size();
	 double diag[size];
	 int nnz = mass_matrix->get_nnz();
	 int row[size]; 
	int col[size+1];
int i;

//#define CHUNKSIZE 1
//#pragma omp parallel shared(diag,row,col) private(i) num_threads(Global<double>::Hermes_omp_get_max_threads())
//{
//#pragma omp for schedule(dynamic, CHUNKSIZE)
	 for(i = 0; i<size; i++){    
	    diag[i] = 0;
	    row[i]= i;
	    col[i]=i;
	 }
//}//Ende parallel

	col[size]=size;// = Anzahl der Nichtnulleintraege (bezieht sich auf theoretisch als naechstes kommenden Eintrag)


	 for(i = 0; i<nnz; i++){    
	    diag[mass_matrix->get_Ai()[i]] += mass_matrix->get_Ax()[i]; 
	 }



	 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
	return lumped_matrix;
}


//artificial Diffusion
UMFPackMatrix<double>* artificialDiffusion_simple(UMFPackMatrix<double>* conv_matrix){
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	double a,b;
	 UMFPackMatrix<double>* diffusion = new UMFPackMatrix<double>;  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0
	for(int i= 0; i<size; i++){
	  for(int j=(i+1);j<size;j++){
	       a= -conv_matrix->get(i,j);
	       b= -conv_matrix->get(j,i);     
	     if((a>=b)&&(a>0.0)){
		 diffusion->add(j,i,a);
		diffusion->add(i,j,a);	
		 diffusion->add(j,j,-a);
		diffusion->add(i,i,-a);	
	     }else if((b>a)&&(b>0.0)){
		 diffusion->add(i,j,b);
		diffusion->add(j,i,b);
	 	 diffusion->add(j,j,-b);
		diffusion->add(i,i,-b);
	     }
	  }
	}
	return diffusion;

}
//Assemble antidiffusive fluxes & Limiter
void antidiffusiveFlux_simple(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix,UMFPackMatrix<double>* conv_matrix,UMFPackMatrix<double>* diffusion,double* u_high, double* u_L, double* u_old, double* flux_double, double* P_plus, double* P_minus, double* Q_plus, double* Q_minus, double* R_plus, double* R_minus  ){
	int ndof = conv_matrix->get_size();
	double alpha,f, plus, minus;
	for(int i=0; i<ndof;i++) flux_double[i]=0.0;
	//Berechnung von P&Q
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
					f = (mass_matrix->get(i,j)+ diffusion->get(i,j)/2.)*(u_high[i]- u_high[j])
									-(mass_matrix->get(i,j)- diffusion->get(i,j)/2.) *(u_old[i]- u_old[j]);	
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
			if(f>0.0)	{ 
				P_plus[i]+=f;
				P_minus[j]-=f;
			}else if (f<0.0){
			 	P_minus[i]+=f;
				P_plus[j]-=f;
			}
			f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i])/time_step; 
			if(f<Q_minus[i]) Q_minus[i] = f;			
			f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j])/time_step;
			if(f>Q_plus[j]) Q_plus[j] = f;	
			if(f<Q_minus[j]) Q_minus[j] = f;						
			
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
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){				
				f = (mass_matrix->get(i,j)+ diffusion->get(i,j)/2.)*(u_high[i]- u_high[j])
					-(mass_matrix->get(i,j)- diffusion->get(i,j)/2.) *(u_old[i]- u_old[j]);	
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
			if(f>0){					
				if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
				else 	alpha = R_plus[i];
			}else{
				if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
				else 	alpha = R_minus[i]; 
			}
			flux_double[i] += alpha*f;
			flux_double[j] -= alpha*f;				
		}
	}
	

}



void lumped_flux_limiter_simple(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix, double* u_L, double* u_H,double* P_plus, double* P_minus, double* Q_plus, double* Q_minus, double* R_plus, double* R_minus){

	int ndof = mass_matrix->get_size();
	UMFPackVector<double>* vec_rhs = new UMFPackVector<double>(ndof);	
		vec_rhs->zero(); 
	double* rhs = new double[ndof];
		lumped_matrix->multiply_with_vector(u_L, rhs); 
	double alpha,f, plus, minus;
	
		//Berechnung von P&Q
		for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
		for(int i =0;i<ndof;i++){
			for(int j =(i+1); j<ndof;j++){		
				f = mass_matrix->get(i,j)*time_step*(u_H[i]- u_H[j]);
			if( (f*(u_H[j]- u_H[i])) > 0.0) f = 0.0; //prelimiting step
				if(f>0.0)	{ 
					P_plus[i]+=f;
					P_minus[j]-=f;
				}else if (f<0.0){
				 	P_minus[i]+=f;
					P_plus[j]-=f;
				}
				f = (u_L[j]-u_L[i]); 
				if(f>Q_plus[i]) Q_plus[i] = f;				
				if(f<Q_minus[i]) Q_minus[i] = f;			
				f= (u_L[i]-u_L[j]); 
				if(f>Q_plus[j]) Q_plus[j] = f;	
				if(f<Q_minus[j]) Q_minus[j] = f;			
			}
		}
		//Berechnung von R
		for(int i=0; i<ndof;i++){
			plus = 1.0; minus = 1.0;		
			if(P_plus[i]!=0.0)  plus = lumped_matrix->get(i,i)*Q_plus[i]/P_plus[i];		
			if(P_minus[i]!=0.0) minus =lumped_matrix->get(i,i)* Q_minus[i]/P_minus[i];			
			if(plus>=1.0) R_plus[i]= 1.0;
			else 	     R_plus[i]= plus;
			if(minus>=1.0) R_minus[i]= 1.0;
			else 	     R_minus[i]= minus;		
		}	
		//Berechnung von alpha & f_i
				alpha = 1.0;
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){				
			f = mass_matrix->get(i,j)*time_step*(u_H[i]- u_H[j]);		
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



