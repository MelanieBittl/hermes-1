int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}



template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* u_prev_time, Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2, CustomWeakFormMassmatrix* massmatrix, HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* al,  std::list<int>* list, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps, int* smooth_elem, double h_start)
{		

	if(!list->empty()) list->clear();


	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	
	int id;

//-----------h-Adapt	
		double* elem_error = new double[space->get_mesh()->get_max_element_id()];
		calc_z_z_error( space, sln, R_h_1, R_h_2, elem_error,smooth_elem);
		double mean_value_z =0.;
		double std_dev_z = 0.0;
		double mean_value_p =0.;
		double std_dev_p = 0.0;
		int counter_h =0; int counter_p=0;
		int counter = space->get_mesh()->get_num_active_elements();

		for_all_active_elements(e, space->get_mesh()){
			if(smooth_elem[e->id]==1){
				mean_value_p += elem_error[e->id];
				std_dev_p += elem_error[e->id]*elem_error[e->id];
				counter_p++;	
			}else{
				mean_value_z += elem_error[e->id];
				std_dev_z += elem_error[e->id]*elem_error[e->id];
				counter_h++;	
			}		
		}	

		std_dev_z -= ((mean_value_z*mean_value_z)/counter_h); mean_value_z/=counter_h;
		std_dev_z/=(counter_h-1);
		std_dev_z = std::sqrt(std_dev_z);

		std_dev_p -= ((mean_value_p*mean_value_p)/counter_p); mean_value_p/=counter_p;
		std_dev_p/=(counter_p-1);
		std_dev_p = std::sqrt(std_dev_p);
		double tol_z = mean_value_z;
		double tol_p = mean_value_p;




		for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if(smooth_elem[e->id]==1){
					if(elem_error[e->id] >tol_p){
										refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;
					}else{
									refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					}
			}else if(elem_error[e->id] >tol_z){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else	if(elem_error[e->id] <EPS){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; i++;
					while((elem_error[e->id]< (EPS/i*1000))&&(i<3)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else {refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++;
			}
		}


		

delete [] elem_error;



Element* elem_neigh=NULL;
	bool p2_neighbor = false;
		for_all_active_elements(e, space->get_mesh()){
			if(elements_to_refine[e->id]==2){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 
								id = elem_neigh->id;	
								if(((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2))||(elements_to_refine[id]==2)){
									p2_neighbor =true;
									break;
								}
							}
					}
					if(p2_neighbor==false){refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++;
			}
					else p2_neighbor = false;			
				}
	}

	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);




//------------------Elemente fuer FCT speichern
	 p2_neighbor =false;
	double elem_diag =0; 
	for_all_active_elements(e, space->get_mesh()){  
		elem_diag=e->get_diameter();
		if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){
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
			if(p2_neighbor==false){list->push_back(e->id);					
			}  // Elemente fuer FCT
			else {p2_neighbor =false;				
			}
		}
	}


	return refine;

}

