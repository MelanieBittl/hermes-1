int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}




	template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,  HPAdapt* adapt, double h_min, double h_max,  int ps)
{		


	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	

				int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   // 0 = nothing..
				int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];	
	memset(elements_to_refine,0, space->get_mesh()->get_max_element_id()*sizeof(int));
memset(no_of_refinement_steps,0, space->get_mesh()->get_max_element_id()*sizeof(int));
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



	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

		delete [] elem_error;
		delete [] elements_to_refine;
		delete [] no_of_refinement_steps;

return refine;


}







