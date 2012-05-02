//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
//Quad hat p-Ordnung fuer horizontal und vertical zusammen space->get_element_order(e->id) = H2D_MAKE_QUAD_ORDER(h-ord, v-ord)
//H2D_GET_H_ORDER(id),  H2D_GET_V_ORDER(id)


template<typename Scalar>
void p1_list_fast(Space<Scalar>* space, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,double* x, double* y){

		Element* e =NULL;
	bool more = false;
	for_all_active_elements(e, space->get_mesh()){
	 //if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){	//Ordnung soll 1 
			space->get_element_assembly_list(e, al);
	  	for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
				int index =  space->get_shapeset()->get_vertex_index(iv);
				Node* vn = e->vn[iv];
				if (space->get_element_order(e->id) == 0) break;
				if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
					for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
						if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
								for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
									if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in liste enhalten
								}
								if(more==false){ dof_list->add_triplet(index, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
																	if(al->get_dof()[j]	!=1.){	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;}
								}
								more = false;
						}
					}
				 }
			}
		//}
	}
}


/*

template<typename Scalar>
void p1_list_fast_neighbor(Space<Scalar>* space, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,std::list<int>* list){

		Element* e =NULL;

	bool more = false;
		std::list<int>::iterator it; it = list->begin(); //Liste sollte nach ids geordnet sein
		int elem_id = *it;
	for_all_active_elements(e, space->get_mesh()){
		 if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){	//Ordnung soll 1 
				if(e->id!=elem_id){
						space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
						int index =  space->get_shapeset()->get_vertex_index(iv);
						Node* vn = e->vn[iv];
						if (space->get_element_order(e->id) == 0) break;
						if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
						for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
									if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
								for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
									if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in liste enhalten
								}
							if(more==false) dof_list->add_triplet(index, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
							more = false;
							}
						}
						 }
					}
				}else{
					if(it!=list->end()){it++; elem_id = *it;}
				}
			}
	}
	for(it=list->begin();it!=list->end();it++){
		elem_id = *it;
	 	if((space->get_element_order(elem_id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(elem_id)==1)){	//Ordnung soll 1 
			e= space->get_mesh()->get_element_fast(elem_id);
			if(e->active){
				space->get_element_assembly_list(e, al);
				for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
						int index =  space->get_shapeset()->get_vertex_index(iv);
						Node* vn = e->vn[iv];
						if (space->get_element_order(elem_id) == 0) break;
						if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
							for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
									if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ //dof=-1 =>dirichlet
										for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
											if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in dof_liste enhalten
										}
										if(more==false){
											for(unsigned int i = 0; i < dof_list_2->get_cnt(); i ++){ 
												if(dof_list_2->get_dof()[i]==al->get_dof()[j]){ dof_list_2->get_idx()[i]=dof_list_2->get_idx()[i]+1;
													more =true; break;}	//ueberpruefen ob dof schon in dof_2_liste enhalten
											}

											 if(more==false){ dof_list_2->add_triplet(1, al->get_dof()[j], 1.0); } //statt e->id, Anzahll an p1 setzen
											
										}
										more = false;
									}
							}
					 	}
					}
			}
		}
	}

}
*/
template<typename Scalar>
void p1_list(Space<Scalar>* space, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al, std::list<int>* list, std::list<int>* list_2,double* x, double* y )
{
	int test =0;
	Element* e =NULL;
	Element* parent =NULL;
	//dof_list->clear();
	//dof_list_2->clear();
	bool more = false;
	int elem_id;
	std::list<int>::iterator it;
  for(it=list->begin();it!=list->end();it++){
			elem_id = *it;
		 	if((space->get_element_order(elem_id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(elem_id)==1)){	//Ordnung soll 1 
				e= space->get_mesh()->get_element_fast(elem_id);
				if(e->active){
					space->get_element_assembly_list(e, al);
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
							int index =  space->get_shapeset()->get_vertex_index(iv);
							Node* vn = e->vn[iv];
							if (space->get_element_order(elem_id) == 0) break;
							if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
								for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
										if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
											for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
												if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in liste enhalten
											}
											if(more==false){ dof_list->add_triplet(e->id, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
																	if(al->get_dof()[j]	!=1.){	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;}
											}
											more = false;
										}
								}
						 	}
						}
				}else{
					parent = e;
					for(int s=0; s<	e->get_num_surf();s++){
							e = parent->sons[s];
							space->get_element_assembly_list(e, al);
							for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
									int index =  space->get_shapeset()->get_vertex_index(iv);
									Node* vn = e->vn[iv];
									if (space->get_element_order(elem_id) == 0) break;
									if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
										for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
												if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
													for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
														if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in liste enhalten
													}
											if(more==false){ dof_list->add_triplet(e->id, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
																	if(al->get_dof()[j]	!=1.){	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;}
											}
													more = false;
												}
										}
								 	}
								}
					}
				}
			}
	}	
	for(it=list_2->begin();it!=list_2->end();it++){
		elem_id = *it;
	 	if((space->get_element_order(elem_id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(elem_id)==1)){	//Ordnung soll 1 
			e= space->get_mesh()->get_element_fast(elem_id);
			if(e->active){
				space->get_element_assembly_list(e, al);
				for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
						int index =  space->get_shapeset()->get_vertex_index(iv);
						Node* vn = e->vn[iv];
						if (space->get_element_order(elem_id) == 0) break;
						if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
							for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
									if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ //dof=-1 =>dirichlet
										for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
											if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in dof_liste enhalten
										}
										if(more==false){
											for(unsigned int i = 0; i < dof_list_2->get_cnt(); i ++){ 
												if(dof_list_2->get_dof()[i]==al->get_dof()[j]){ dof_list_2->get_idx()[i]=dof_list_2->get_idx()[i]+1;
													more =true; break;}	//ueberpruefen ob dof schon in dof_2_liste enhalten
											}

	 											if(more==false){ dof_list->add_triplet(e->id, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
																	if(al->get_dof()[j]	!=1.){	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;}
											}
											
										}
										more = false;
									}
							}
					 	}
					}
			}else{
					parent = e;
					for(int s=0; s<	e->get_num_surf();s++){
						e = parent->sons[s];
						space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
								int index =  space->get_shapeset()->get_vertex_index(iv);
								Node* vn = e->vn[iv];
								if (space->get_element_order(elem_id) == 0) break;
								if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
									for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
											if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ //dof=-1 =>dirichlet
												for(unsigned int i = 0; i < dof_list->get_cnt(); i ++){ 
													if(dof_list->get_dof()[i]==al->get_dof()[j]){ more =true; break;}	//ueberpruefen ob dof schon in dof_liste enhalten
												}
												if(more==false){
													for(unsigned int i = 0; i < dof_list_2->get_cnt(); i ++){ 
														if(dof_list_2->get_dof()[i]==al->get_dof()[j]){ dof_list_2->get_idx()[i]=dof_list_2->get_idx()[i]+1;
															more =true; break;}	//ueberpruefen ob dof schon in dof_2_liste enhalten
													}

													// if(more==false){ dof_list_2->add_triplet(1, al->get_dof()[j], 1.0);  test++;} //statt e->id, Anzahll an p1 setzen
												if(more==false){ dof_list->add_triplet(e->id, al->get_dof()[j], 1.0);  //dof=-1 =>dirichlet
																	if(al->get_dof()[j]	!=1.){	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;}
											}
																							//=>haengender Knoten zwischen Problemknoten!
												}
												more = false;
											}
									}
							 	}
							}
						}
			}
		}
	}		



}


