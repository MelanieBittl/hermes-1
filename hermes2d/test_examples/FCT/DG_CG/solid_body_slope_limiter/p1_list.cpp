
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
	if(elem_diag>h_start){	p2_neighbor =true;  //vergroebert!
				}else{
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
					{  
					  // Class for finding all neighbors of the element e on the mesh space->get_mesh().
					  NeighborSearch<double> ns(e, space->get_mesh());
					  // Set this to ignore errors of the edge iv being a boundary edge (which has no neighbors).
					  ns.set_ignore_errors(true);
					  // Look for the neighbors over the edge iv.
					  ns.set_active_edge(iv);
					  // Iterate through the found neighbors.
					  for(unsigned int i = 0; i < ns.get_neighbors()->size(); i++)
					  {
						 id = ns.get_neighbors()->at(i)->id;
						 if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2))
						 {
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


