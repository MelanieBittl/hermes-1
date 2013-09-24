//artificial Diffusion
CSCMatrix<double>* artificialDiffusion(CSCMatrix<double>* conv_matrix){
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	double a,b;
	 CSCMatrix<double>* diffusion = new CSCMatrix<double>;  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0

	double* Ax = conv_matrix->get_Ax();
	int* Ai = conv_matrix->get_Ai();
	int* Ap = conv_matrix->get_Ap();

		for(int j = 0; j<size; j++){ //Spalten durchlaufen
				for(int indx = Ap[j]; indx<Ap[j+1];indx++){	
							int i = Ai[indx];	
							if(j<i){
								a= -Ax[indx];
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
	}


	/*for(int i= 0; i<size; i++){
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
	}*/
	return diffusion;

}



  //artificial Diffusion
CSCMatrix<double>* artificialDiffusion(bool* fct, CSCMatrix<double>* conv_matrix){
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	double a,b;
	 CSCMatrix<double>* diffusion = new CSCMatrix<double>;  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0

	double* Ax = conv_matrix->get_Ax();
	int* Ai = conv_matrix->get_Ai();
	int* Ap = conv_matrix->get_Ap();

		for(int j = 0; j<size; j++){ //Spalten durchlaufen
				if(fct[j]== false) continue;
				for(int indx = Ap[j]; indx<Ap[j+1];indx++){	
							int i = Ai[indx];	
							if(fct[i]== false) continue;
							if(j<i){
								a= -Ax[indx];
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
	}


	return diffusion;

}

//artificial Diffusion
template<typename Scalar>
CSCMatrix<Scalar>* artificialDiffusion(AsmList<Scalar>* al,CSCMatrix<Scalar>* conv_matrix)
{ //al=NULL => diffusion=0
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	Scalar a,b;
	 CSCMatrix<Scalar>* diffusion = new CSCMatrix<Scalar>;  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0
	if(al!=NULL){ //nur wenn Liste fuer Vertex-DOFS vorhanden
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
		 for(unsigned int j = (i+1); j < al->get_cnt(); j ++){
			if(al->get_dof()[i]!=al->get_dof()[j]){		
						     a= -conv_matrix->get(al->get_dof()[i],al->get_dof()[j]);
						    	b= -conv_matrix->get(al->get_dof()[j],al->get_dof()[i]);     
						   if((a>=b)&&(a>0.0)){
						 diffusion->add(al->get_dof()[j],al->get_dof()[i],a);
						diffusion->add(al->get_dof()[i],al->get_dof()[j],a);	
						 diffusion->add(al->get_dof()[j],al->get_dof()[j],-a);
						diffusion->add(al->get_dof()[i],al->get_dof()[i],-a);	
						   }else if((b>a)&&(b>0.0)){
						 diffusion->add(al->get_dof()[i],al->get_dof()[j],b);
						diffusion->add(al->get_dof()[j],al->get_dof()[i],b);
					 	 diffusion->add(al->get_dof()[j],al->get_dof()[j],-b);
						diffusion->add(al->get_dof()[i],al->get_dof()[i],-b);
			   				}
				}
			}
		}
	}



	return diffusion;

}
