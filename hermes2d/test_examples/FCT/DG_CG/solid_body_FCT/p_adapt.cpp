int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}



template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2, CustomWeakFormMassmatrix* massmatrix, HPAdapt* adapt, AsmList<Scalar>* al,  double h_min, double h_max, int ts, int ps, int* smooth_elem, double h_start)
{		

	int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   // 0 = nothing..
	int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];	
	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	
	int id;

//-----------h-Adapt	
		double* elem_error = new double[space->get_mesh()->get_max_element_id()];
		calc_z_z_error( space, sln, R_h_1, R_h_2, elem_error);
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



	for_all_active_elements(e, space->get_mesh())
	{
		no_of_refinement_steps[e->id]=0;	elements_to_refine[e->id] = 0;
			int i = 1;	
		/*	if(elem_error[e->id] <EPS){ //vergroebern
							refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]< (EPS/i*100))&&(i<4)){
					 	no_of_refinement_steps[e->id]++; i++;
					}
			}else */
			if(smooth_elem[e->id]==1){  //glatt => p erhoehen
				//	if(elem_error[e->id] >tol_p){
										refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;
										/*	while((elem_error[e->id]> (tol_p+i*std_dev_p))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
										} */
					//}
					/*else if(elem_error[e->id] >tol_z){ // h verkleinern
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
						while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
						 													no_of_refinement_steps[e->id]++; i++;
						}
					}*/
		/*	}else if(elem_error[e->id] >tol_z){ //nichtglatt => h verkleinern
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}*/
			}
		}



delete [] elem_error;


	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);



			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
	return refine;

}

