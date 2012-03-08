

//Mass lumping an den Stellen von AsmList, sonst standard Massmatrix
template<typename Scalar>
UMFPackMatrix<Scalar>* massLumping(AsmList<Scalar>* al,AsmList<Scalar>* al_2, UMFPackMatrix<Scalar>* mass_matrix)
{  //al=NULL=>lumped=mass
	 UMFPackMatrix<Scalar>* lumped_matrix = new UMFPackMatrix<Scalar>;   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	Scalar a =0.0;
	if(al!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
			for(unsigned int j = (i+1); j < al->get_cnt(); j ++){	//Mc_ij = Mc_ji				
				if(al->get_dof()[i]!=al->get_dof()[j]){  //kein Diagonaleintrag!
					if((al->get_dof()[i]>=size)||(al->get_dof()[j]>=size)) error("massLumping:wrong DOF");				
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
	if(al_2!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
		for(unsigned int i = 0; i < al->get_cnt(); i ++){
			for(unsigned int j = 0; j < al_2->get_cnt(); j ++){	//Mc_ij = Mc_ji		
				if(al->get_dof()[i]!=al_2->get_dof()[j]){  //kein Diagonaleintrag!
					if((al->get_dof()[i]>=size)||(al_2->get_dof()[j]>=size)) error("massLumping:wrong DOF");				
					a = lumped_matrix->get(al->get_dof()[i],al_2->get_dof()[j]);
					if(a!=0.0){
						lumped_matrix->add(al->get_dof()[i],al->get_dof()[i],a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(al_2->get_dof()[j],al_2->get_dof()[j],a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(al->get_dof()[i],al_2->get_dof()[j],-a);	//i,j Eintrag auf 0 setzen
						lumped_matrix->add(al_2->get_dof()[j],al->get_dof()[i],-a);	//i,j Eintrag auf 0 setzen
					}				
				}
			}			
		}
	}
	
	return lumped_matrix;
}

//Mass lumping an den Stellen von AsmList, sonst 0
template<typename Scalar>
UMFPackMatrix<Scalar>* pure_p1_massLumping(AsmList<Scalar>* al,UMFPackMatrix<Scalar>* mass_matrix)
{  //al=NULL=>lumped=mass
	 UMFPackMatrix<Scalar>* lumped_matrix = new UMFPackMatrix<Scalar>;   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	lumped_matrix->zero();
	Scalar a =0.0;
	if(al!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
	for(unsigned int i = 0; i < al->get_cnt(); i ++){
		for(unsigned int j = i; j < al->get_cnt(); j ++){	//Mc_ij = Mc_ji				
				if((al->get_dof()[i]>=size)||(al->get_dof()[j]>=size)) error("massLumping:wrong DOF");				
				a = mass_matrix->get(al->get_dof()[i],al->get_dof()[j]);
				if(a!=0.0){
					lumped_matrix->add(al->get_dof()[i],al->get_dof()[i],a);    //zur Diagonale hinzufuegen
					if(al->get_dof()[i]!=al->get_dof()[j]) lumped_matrix->add(al->get_dof()[j],al->get_dof()[j],a);    //zur Diagonale hinzufuegen
				}				
			
		}			
	}
	}
	
	return lumped_matrix;
}
