//Mass lumping
CSCMatrix<double>* massLumping(CSCMatrix<double>* mass_matrix){
	 CSCMatrix<double>* lumped_matrix = new CSCMatrix<double>;   //M_L
	 int size = mass_matrix->get_size();
	 double diag[size];
	 int nnz = mass_matrix->get_nnz();
	 int row[size]; 
	int col[size+1];
int i;
	 for(i = 0; i<size; i++){    
	    diag[i] = 0;
	    row[i]= i;
	    col[i]=i;
	 }
	col[size]=size;// = Anzahl der Nichtnulleintraege (bezieht sich auf theoretisch als naechstes kommenden Eintrag)
	 for(i = 0; i<nnz; i++){    
	    diag[mass_matrix->get_Ai()[i]] += mass_matrix->get_Ax()[i]; 
	 }
	 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
	return lumped_matrix;
}

//Mass lumping an den Stellen von AsmList, sonst standard Massmatrix
template<typename Scalar>
CSCMatrix<Scalar>* massLumping(bool* fct, CSCMatrix<Scalar>* mass_matrix)
{  //al=NULL=>lumped=mass
	 CSCMatrix<Scalar>* lumped_matrix = new CSCMatrix<Scalar>;   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	double* Ax = lumped_matrix->get_Ax();
	int* Ai = lumped_matrix->get_Ai();
	int* Ap = lumped_matrix->get_Ap();
	Scalar a =0.0;

		for(int j = 0; j<size; j++){ //Spalten durchlaufen
				if(fct[j]== false) continue;
				for(int indx = Ap[j]; indx<Ap[j+1];indx++){	
							int i = Ai[indx];	
							if(fct[i]== false) continue;
							if(j<i){
								a = Ax[indx];
								if(a!=0.0){
									lumped_matrix->add(i,i,a);    //zur Diagonale hinzufuegen
									lumped_matrix->add(j,j,a);    //zur Diagonale hinzufuegen
									lumped_matrix->add(i,j,-a);	  //i,j Eintrag auf 0 setzen
									lumped_matrix->add(j,i,-a);	 //j,i Eintrag auf 0 setzen
								}	
						}
				}
			}
	return lumped_matrix;
}



//Mass lumping an den Stellen von AsmList, sonst standard Massmatrix
template<typename Scalar>
CSCMatrix<Scalar>* massLumping(AsmList<Scalar>* al, CSCMatrix<Scalar>* mass_matrix)
{  //al=NULL=>lumped=mass
	 CSCMatrix<Scalar>* lumped_matrix = new CSCMatrix<Scalar>;   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	Scalar a =0.0;
	if(al!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
			for(unsigned int j = (i+1); j < al->get_cnt(); j ++){	//Mc_ij = Mc_ji				
				if(al->get_dof()[i]!=al->get_dof()[j]){  //kein Diagonaleintrag!
					if((al->get_dof()[i]>=size)||(al->get_dof()[j]>=size)) printf("massLumping:wrong DOF");				
					a = lumped_matrix->get(al->get_dof()[i],al->get_dof()[j]);
					if(a!=0.0){
						lumped_matrix->add(al->get_dof()[i],al->get_dof()[i],a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(al->get_dof()[j],al->get_dof()[j],a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(al->get_dof()[i],al->get_dof()[j],-a);	//i,j Eintrag auf 0 setzen
						lumped_matrix->add(al->get_dof()[j],al->get_dof()[i],-a);	//i,j Eintrag auf 0 setzen
					}				
				}
			}			
		}
	}	
	return lumped_matrix;
}

