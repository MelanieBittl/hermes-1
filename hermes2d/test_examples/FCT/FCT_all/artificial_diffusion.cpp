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


	return diffusion;

}

