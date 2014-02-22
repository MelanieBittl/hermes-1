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
	    diag[i] = 0.;
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

//Mass lumping an den Stellen von fct, sonst standard Massmatrix
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

bool* get_vertex_dofs(Hermes::vector<SpaceSharedPtr<double> > spaces)
{

	SpaceSharedPtr<double> space = spaces[0];
	int ndof = Space<double>::get_num_dofs(spaces);
	int dof_rho = space->get_num_dofs();
	bool* fct= new bool[ndof];
	for(int i =0; i<ndof;i++) fct[i] = false;
	Element* e;	AsmList<double> al;
	for_all_active_elements(e, space->get_mesh())
		{

			space->get_element_assembly_list(e, &al);
			for(unsigned int j = 0; j < al.get_cnt(); j ++)
			{			 
			  	for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
				  int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
						Node* vn = e->vn[iv];
						if (space->get_element_order(e->id) == 0) break;
						if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
							for(unsigned int j = 0; j < al.get_cnt(); j ++){			 
								if((al.get_idx()[j]==index)&&(al.get_dof()[j]!=-1.0)){ 
												int dof = al.get_dof()[j];
												for(int k = 0; k<4; k++)
														fct[dof+k*dof_rho]=true;

								}
							}
						 }
					}
			}




		
				
								
		}
 return fct;
}


