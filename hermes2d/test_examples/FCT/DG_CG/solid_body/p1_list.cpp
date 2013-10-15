//Werte von u_prev_time in dofs speichern
template<typename Scalar>
void coord_dof(SpaceSharedPtr<Scalar> space, AsmList<Scalar>* al,double* x, MeshFunctionSharedPtr<Scalar>  sln){

		Element* e =NULL;
	bool more = false;
	for_all_active_elements(e, space->get_mesh()){
	 //if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){	//Ordnung soll 1 
			space->get_element_assembly_list(e, al);
	  	for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		  int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
				Node* vn = e->vn[iv];
				if (space->get_element_order(e->id) == 0) break;
				if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
					for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
						if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
										if(x[al->get_dof()[j]]==0.)		x[al->get_dof()[j]]	= sln->get_pt_value(vn->x, vn->y);
						}
					}
				 }
			}
		//}
	}
}

template<typename Scalar>
void p1_list(SpaceSharedPtr<Scalar> space, bool* fct, AsmList<Scalar>* al, double* x, MeshFunctionSharedPtr<Scalar>  sln, double h_start )
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
				if(elem_diag>h_start){	p2_neighbor =true;
				}else{
							for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
								 	elem_neigh = e->get_neighbor(iv);
									if(elem_neigh!=NULL){ 
											 id = elem_neigh->id;	
											if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
												p2_neighbor =true; 
												break;
											}else if(elem_diag != elem_neigh->get_diameter()){
														// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
												p2_neighbor =true; 
												break;
											}
									}
									if(e->vn[iv]->is_constrained_vertex() ==true)	{	p2_neighbor =true; 
													break;
									}
							}
				}
				if(p2_neighbor==false){
							space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		  int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
								Node* vn = e->vn[iv];
								if (space->get_element_order(elem_id) == 0) break;
								if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
									for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
											if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
													if(fct[al->get_dof()[j]]==false){
																			x[al->get_dof()[j]]	= sln->get_pt_value(vn->x, vn->y);
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

template<typename Scalar>
void p1_list(SpaceSharedPtr<Scalar> space, bool* fct, AsmList<Scalar>* al, double h_start )
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
				if(elem_diag>h_start){	p2_neighbor =true;
				}else{
							for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
								 	elem_neigh = e->get_neighbor(iv);
									if(elem_neigh!=NULL){ 
											 id = elem_neigh->id;	
											if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
												p2_neighbor =true; 
												break;
											}else if(elem_diag != elem_neigh->get_diameter()){
														// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
												p2_neighbor =true; 
												break;
											}
									}
									if(e->vn[iv]->is_constrained_vertex() ==true)	{	p2_neighbor =true; 
													break;
									}
							}
				}
				if(p2_neighbor==false){
							space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		  int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
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

template<typename Scalar>
void vertex_dof_list(SpaceSharedPtr<Scalar> space,bool* fct,  AsmList<Scalar>* al )
{
	Element* e =NULL;
	Element* elem_neigh=NULL;
	bool more = false;
	bool  p2_neighbor =false;
	double elem_diag =0; 
	int elem_id,id;
	for_all_active_elements(e, space->get_mesh()){ 
				elem_diag=e->get_diameter();
				elem_id= e->id;
				space->get_element_assembly_list(e, al);
				for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		 			 int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_TRIANGLE);
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
	} 
		
}





