int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}


	template<typename Scalar>
bool init_coarse(Space<Scalar>* space,Solution<Scalar>* sln, HPAdapt* adapt, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max){

	bool refine = false;  //->verfeinern oder vergroebern?	
	Element* e = NULL;
	int dof = space->get_num_dofs();
	Scalar* zeros = new Scalar[dof]; for(int i =0; i<dof; i++) zeros[i] = 0.;
		Solution<Scalar>* zero_sln = new Solution<Scalar>(space, zeros);
			Adapt<Scalar>* adaptivity_h = new Adapt<Scalar>(space ,HERMES_H1_NORM);
		double err_est_h = adaptivity_h->calc_err_est(zero_sln,sln,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
		for_all_active_elements(e, space->get_mesh()){
		 if((std::sqrt(adaptivity_h->get_element_error_squared(0, e->id))<EPS)){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++;
			}
		}
	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);


	return refine;
}


	template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,  std::list<int>* list,std::list<int>* neighbor, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps)
//bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,  HPAdapt* adapt, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max,  int ps)
{		


	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	


		double* elem_error = new double[space->get_mesh()->get_max_element_id()];

		calc_z_z_error( space, sln, R_h_1, R_h_2, elem_error);

		double mean_value_h=0.;
		double std_dev_h = 0.0;
		int counter = space->get_mesh()->get_num_active_elements();	
		for_all_active_elements(e, space->get_mesh()){
				mean_value_h += elem_error[e->id];
				std_dev_h += elem_error[e->id]*elem_error[e->id];			
		}	

		std_dev_h -= ((mean_value_h*mean_value_h)/counter); mean_value_h/=counter;
		std_dev_h/=(counter-1);
		std_dev_h = std::sqrt(std_dev_h);
		double tol_h = mean_value_h;


		for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if(elem_error[e->id] >tol_h){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_h+i*std_dev_h))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else	if(elem_error[e->id] <EPS){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]< (EPS/i*100))&&(i<3)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}
		}	

		delete [] elem_error;

	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);


//------------------Elemente fuer FCT speichern
/*	if(!list->empty()) list->clear();
	if(!neighbor->empty()) neighbor->clear();
	bool p2_neighbor =false;Element* elem_neigh =NULL;
	double elem_diag =0; 
	int id;
	for_all_active_elements(e, space->get_mesh()){  
		elem_diag=e->get_diameter();
			for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 
								id = elem_neigh->id;	
							 if(elem_diag != elem_neigh->get_diameter()){
											// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
									p2_neighbor =true; //neighbor->push_back(e->id);
									break;
								}
						}
						if(e->vn[iv]->is_constrained_vertex() ==true)	{	p2_neighbor =true; //neighbor->push_back(e->id);
										break;
						}
			}
			if(p2_neighbor==false){list->push_back(e->id);					
			}  // Elemente fuer FCT
			else {p2_neighbor =false;				
			}
		
	}*/

return refine;


}







