#include "hp_adapt.h"

bool refine_elem(Hermes::vector<Space<double>*> spaces, Element* e, int ref, double h_min){
			bool refined = true;
			int order = spaces[0]->get_element_order(e->id);
			if(e->get_diameter()>h_min){
				if (e->active)
		  			spaces[0]->get_mesh()->refine_element_id(e->id);				
				for (int j = 0; j < 4; j++){
							for(unsigned int i = 0; i < spaces.size(); i++)
									spaces[i]->set_element_order_internal(e->sons[j]->id, order);	

				}	
				if(ref>1) {
					for (int j = 0; j < 4; j++)
						refined = refine_elem(spaces, e->sons[j], ref-1, h_min);
				}
			}		
		return refined;
}



bool refine_elem_order(Hermes::vector<Space<double>*> spaces, Element* e, int ref, int order, double h_min){
			bool refined = true;
			if(e->get_diameter()>h_min){			
				if (e->active)
		  			spaces[0]->get_mesh()->refine_element_id(e->id);				
				for (int j = 0; j < 4; j++){
						for(unsigned int i = 0; i < spaces.size(); i++)
									spaces[i]->set_element_order_internal(e->sons[j]->id, order);	
				}	
				if(ref>1) {
					for (int j = 0; j < 4; j++)
						refined = refine_elem(spaces, e->sons[j], ref-1, h_min);
				}
			}
		return refined;
}


bool coarse_elem(Hermes::vector<Space<double>*> spaces, Element* e, int ref, int order, double h_max, int* elements_to_refine){
				bool coarsed =true;
				Element* parent = e->parent;				
				if(parent!=NULL) {	
					if(parent->get_diameter()< h_max){
						bool active_sons = true;	
						for (int j = 0; j < 4; j++){ 
							if(parent->sons[j]!=NULL){
								if(!parent->sons[j]->active){ active_sons = false; break;
								}else{
									if(elements_to_refine[parent->sons[j]->id]==4)elements_to_refine[parent->sons[j]->id] =0;
									if(elements_to_refine[parent->sons[j]->id]==1){active_sons = false; break;}
								}
							}
						}
										
						if(active_sons==true){		
								for(unsigned int i = 0; i < spaces.size(); i++)
									spaces[i]->set_element_order_internal(parent->id, order);						
							//space->set_element_order_internal(parent->id,order);
							spaces[0]->get_mesh()->unrefine_element_id(parent->id);								
						}
						if(ref>1){
							coarsed = coarse_elem(spaces, parent, ref-1,order, h_max, elements_to_refine);
						}
					}
				}
			return coarsed;
}




bool HPAdapt::adapt(int* elements_to_refine,int* no_of_refinement_steps, int* smooth_elem, int regularize)
{ 
	if(elements_to_refine==NULL) return false;

	bool changed = false;
	
	Space<double>* space = this->spaces[0];
	int order, v_ord, h_ord,ref;
	int n_dof = space->get_num_dofs();
	//apply refinements
	Element* e;	
	for_all_active_elements(e, space->get_mesh()){	
		if((elements_to_refine[e->id] == 1)&&(n_dof<max_dof)){				//h refine 
			if((smooth_elem!=NULL)&&(smooth_elem[e->id] == 1)) //smooth => p increase	
			{						
				if(e->is_triangle()==true)
				{
					order = space->get_element_order(e->id); 
					order++;
					if(order>max_p) order = max_p;
					space->set_element_order_internal(e->id, order);
				}
				else if(e->is_quad()==true)
				{
					v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
					h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
					h_ord++; if(h_ord>max_p) h_ord =max_p;
					v_ord++; if(v_ord>max_p) v_ord =max_p;
					order = H2D_MAKE_QUAD_ORDER(h_ord, v_ord);
				for(unsigned int i = 0; i < spaces.size(); i++)
					spaces[i]->set_element_order_internal(e->id, order);
					changed = true;
				}					
			}else{		
				ref = no_of_refinement_steps[e->id];
				changed = refine_elem(spaces, e, ref,h_min);
				}
		}else if(elements_to_refine[e->id] == 4){				//h coarse 	
				ref = no_of_refinement_steps[e->id];
				//order = space->get_element_order(e->id);	
				order = H2D_MAKE_QUAD_ORDER(1, 1);	
				changed= coarse_elem(spaces, e, ref,order, h_max, elements_to_refine);					
		}
	}

	  
	if(changed==false){ info("nothing to refine/coarse");return false;}


   // since space changed, assign dofs:
      for(unsigned int i = 0; i < this->spaces.size(); i++)
					this->spaces[i]->assign_dofs();



  return true;
}








