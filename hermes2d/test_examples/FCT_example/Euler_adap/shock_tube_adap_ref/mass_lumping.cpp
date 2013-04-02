//Mass lumping
UMFPackMatrix<double>* massLumping(UMFPackMatrix<double>* mass_matrix){
	 UMFPackMatrix<double>* lumped_matrix = new UMFPackMatrix<double>;   //M_L
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

//Mass lumping an den Stellen von fct, sonst standard Massmatrix
template<typename Scalar>
UMFPackMatrix<Scalar>* massLumping(bool* fct, UMFPackMatrix<Scalar>* mass_matrix)
{  //al=NULL=>lumped=mass
	 UMFPackMatrix<Scalar>* lumped_matrix = new UMFPackMatrix<Scalar>;   //M_L
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



template<typename Scalar>
void p1_list(Space<Scalar>* space, bool* fct, AsmList<Scalar>* al)
{
	Element* e =NULL;
Element* elem_neigh=NULL;
	bool more = false;
	bool  p2_neighbor =false;
	double elem_diag =0; 
	int elem_id,id;
	for_all_active_elements(e, space->get_mesh()){  
		if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){
				elem_diag=e->get_diameter();
				elem_id= e->id;
				//if(elem_diag>h_start){	p2_neighbor =true;
				//}else{
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 
							 id = elem_neigh->id;	
							if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
								p2_neighbor =true; 
								break;
							}else if(elem_diag != elem_neigh->get_diameter()){//Nachbar ist kleiner=> haengender Knoten, kein FCT ueber haengendne Knoten 
								p2_neighbor =true; 
								break;
							}
						}
						if(e->vn[iv]->is_constrained_vertex() ==true)	{	p2_neighbor =true; 
										break;
						}
					}
				//}
				if(p2_neighbor==false){
							space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
								int index =  space->get_shapeset()->get_vertex_index(iv);
								Node* vn = e->vn[iv];
								if (space->get_element_order(elem_id) == 0) break;
								if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
									for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
											if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
													if(fct[al->get_dof()[j]]==false){
																		fct[al->get_dof()[j]]=true;
													}
											}
									}
							 	}
							}					
				}  // Elemente fuer FCT
					else {p2_neighbor =false;				
				}
		}
	}

}


