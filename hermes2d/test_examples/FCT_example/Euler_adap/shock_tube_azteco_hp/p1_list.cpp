
void p1_list(Space<double>* space, bool* fct, AsmList<double>* al,int dof_rho, int dof_v_x, int dof_v_y,int dof_e)
{
int ndof = space->get_num_dofs();
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
				for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
				 	elem_neigh = e->get_neighbor(iv);
					if(elem_neigh!=NULL)
					{ 
						 id = elem_neigh->id;	
						if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
							p2_neighbor =true; 
							break;
						}else if(elem_diag != elem_neigh->get_diameter()){//Nachbar ist kleiner=> haengender Knoten, kein FCT ueber haengendne Knoten 
							p2_neighbor =true; 
							break;
						}
					}
					if(e->vn[iv]->is_constrained_vertex() ==true)	
					{			
							p2_neighbor =true; 
							break;
					}
				}

				if(p2_neighbor==false){
							space->get_element_assembly_list(e, al);
						for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
								int index =  space->get_shapeset()->get_vertex_index(iv);
								Node* vn = e->vn[iv];
								if (space->get_element_order(elem_id) == 0) break;
								if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
									for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
											if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0)){ 
													if(al->get_dof()[j]>dof_rho) error("p1_List!\n");
													if(fct[al->get_dof()[j]]==false){
																		fct[al->get_dof()[j]]=true;
																		fct[al->get_dof()[j]+dof_rho]=true;
																		fct[al->get_dof()[j]+dof_rho+dof_v_x]=true;
																		fct[al->get_dof()[j]+dof_rho+dof_v_x+dof_v_y]=true;
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


