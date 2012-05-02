//artificial Diffusion
UMFPackMatrix<double>* artificialDiffusion(UMFPackMatrix<double>* conv_matrix){
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

//artificial Diffusion
template<typename Scalar>
UMFPackMatrix<Scalar>* artificialDiffusion(AsmList<Scalar>* al,UMFPackMatrix<Scalar>* conv_matrix)
{ //al=NULL => diffusion=0
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	Scalar a,b;
	 UMFPackMatrix<Scalar>* diffusion = new UMFPackMatrix<Scalar>;  
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
