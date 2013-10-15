#include "hp_adapt.h"

bool refine_elem(SpaceSharedPtr<double> space, Element* e, int ref, double h_min){
			bool refined = true;
			int order = space->get_element_order(e->id);
			if(e->get_diameter()>h_min){
				if (e->active)
		  			space->get_mesh()->refine_element_id(e->id);				
				for (int j = 0; j < 4; j++){
				  space->set_element_order_internal(e->sons[j]->id, order);
				}	
				if(ref>1) {
					for (int j = 0; j < 4; j++)
						refined = refine_elem(space, e->sons[j], ref-1, h_min);
				}
			}		
		return refined;
}



bool coarse_elem(SpaceSharedPtr<double> space, Element* e, int ref, int order, double h_max, int* elements_to_refine){
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
							space->set_element_order_internal(parent->id,order);
							space->get_mesh()->unrefine_element_id(parent->id);								
						}
						if(ref>1){
							coarsed = coarse_elem(space, parent, ref-1,order, h_max, elements_to_refine);
						}
					}
				}
			return coarsed;
}


bool HPAdapt::adapt(int* elements_to_refine,int* no_of_refinement_steps, int max_p,double h_min,double h_max,int max_dof, int regularize)
{ 
	if(elements_to_refine==NULL) return false;
	  //get meshes
  MeshSharedPtr meshes[H2D_MAX_COMPONENTS];
  for (int j = 0; j < this->num; j++) {
    meshes[j] = this->spaces[j]->get_mesh();
  }
	bool changed = false;
	SpaceSharedPtr<double> space = this->spaces[0];
	int order, v_ord, h_ord,ref;
	int n_dof = space->get_num_dofs();
  //apply refinements
	Element* e;	
	for_all_active_elements(e, space->get_mesh()){
		if((elements_to_refine[e->id] == 1)&&(n_dof<max_dof)){				//h refine 
				ref = no_of_refinement_steps[e->id];
				changed = refine_elem(space, e, ref,h_min);
		}else if(elements_to_refine[e->id] == 4){				//h coarse 	
				ref = no_of_refinement_steps[e->id];
				order = space->get_element_order(e->id);		
				changed= coarse_elem(space, e, ref,order, h_max, elements_to_refine);			
		
		}		
	}
	  
	if(changed==false){ info("nothing to refine/coarsen");return false;}



   // since space changed, assign dofs:
      for(unsigned int i = 0; i < this->spaces.size(); i++)
        this->spaces[i]->assign_dofs();

  return true;
}








