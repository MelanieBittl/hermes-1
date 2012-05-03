//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
//Quad hat p-Ordnung fuer horizontal und vertical zusammen space->get_element_order(e->id) = H2D_MAKE_QUAD_ORDER(h-ord, v-ord)
//H2D_GET_H_ORDER(id),  H2D_GET_V_ORDER(id)


template<typename Scalar>
void p1_list_fast(Space<Scalar>* space, AsmList<Scalar>* dof_list, AsmList<Scalar>* al,double* x, double* y){
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
																	x[al->get_dof()[j]]	= vn->x;y[al->get_dof()[j]]	= vn->y;
								}
								more = false;
						}
					}
				 }
			}
		//}
	}
}

