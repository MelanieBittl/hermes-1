int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}



template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2, HPAdapt* adapt, AsmList<Scalar>* al,  double h_min, double h_max, int ts, int ps, int* smooth_elem, double h_start)
{		

	int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   // 0 = nothing..
	int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];	
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
/*
for_all_active_elements(e, space->get_mesh()){
		no_of_refinement_steps[e->id]=0;	elements_to_refine[e->id] = 0;
		if(smooth_elem[e->id]==1){refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;}
		}*/

	for_all_active_elements(e, space->get_mesh()){
		no_of_refinement_steps[e->id]=0;	elements_to_refine[e->id] = 0;
			int i = 1;	
			if(elem_error[e->id] <EPS_h){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]< (EPS_h/i*100))&&(i<4)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else if(smooth_elem[e->id]==1){
					if(elem_error[e->id] >tol_p){
										refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;
						while((elem_error[e->id]> (tol_p+i*std_dev_p))&&(i<(P_MAX-1))){
						 													no_of_refinement_steps[e->id]++; i++;
						}
				}
				/*/else if(elem_error[e->id] >EPS_h){
				else if(elem_error[e->id] >tol_z){
									refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
														while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
														}
					}*/
			}else if(elem_error[e->id] >tol_z){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else{
				//refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
			}
		}
		

delete [] elem_error;


/*
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
								//	break;
								}
									else{
											refine = true; elements_to_refine[id] = 1; no_of_refinement_steps[id]=1;
								}
							}
					}
					if(p2_neighbor==false){refine = true; elements_to_refine[e->id] = 1; 
					}	else p2_neighbor = false;			
				}
	}
*/

Element* elem_neigh=NULL;Element* elem_neigh_2=NULL; int id_2;
	bool p2_neighbor = false;
		for_all_active_elements(e, space->get_mesh()){
			if(elements_to_refine[e->id]==2){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 
								id = elem_neigh->id;	
								if(((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2))||(elements_to_refine[id]==2)){
									p2_neighbor =true;
									//break;
								}else{
										elements_to_refine[id] = 1; no_of_refinement_steps[id]=1; 
										for (unsigned int ivn = 0; ivn < elem_neigh->get_nvert(); ivn++){  
											 	elem_neigh_2 = elem_neigh->get_neighbor(ivn);
												if(elem_neigh_2!=NULL){ 
														id_2 = elem_neigh_2->id;	
														if(elements_to_refine[id_2]==0){
																elements_to_refine[id_2] = 1; no_of_refinement_steps[id_2]=1; 
														}
													}
											}



								}
							}
					}
					if(p2_neighbor==false){refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[id]=1; 
			}
					else p2_neighbor = false;			
				}
	}

	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);



			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
	return refine;

}

