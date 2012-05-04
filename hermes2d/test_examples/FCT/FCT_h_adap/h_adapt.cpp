int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}


	template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,  HPAdapt* adapt,double h_min, double h_max)
{		

	int ndof = space->get_num_dofs();
	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	
	int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   // 0 = nothing..
	int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];	
		double* elem_error = new double[space->get_mesh()->get_max_element_id()];

//-------------- aus reg_estimator------------------------------
	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(sln);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(sln);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ; 
	UMFPackVector<double> * rhs_1 = new UMFPackVector<double>(ndof);
	UMFPackVector<double> * rhs_2 = new UMFPackVector<double>(ndof);
	dp_1->assemble(mass_matrix,rhs_1); 
	dp_2->assemble(rhs_2);
	UMFPackLinearSolver<double> * solver_1 = new UMFPackLinearSolver<double> (mass_matrix,rhs_1);
	if(solver_1->solve()){ 				
		Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , space, R_h_1);	
	}else error ("Matrix solver failed.\n");
	UMFPackLinearSolver<double> * solver_2 = new UMFPackLinearSolver<double> (mass_matrix,rhs_2);	
	if(solver_2->solve()){ 				
		Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , space, R_h_2);	
	}else error ("Matrix solver failed.\n");
//--------------------------


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
			no_of_refinement_steps[e->id]=0;					
			int i = 1;	
			if(elem_error[e->id] >tol_h){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_h+i*std_dev_h))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else	if(elem_error[e->id] <EPS){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]< (EPS/i*100))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else {
					elements_to_refine[e->id] = 0;
			}
		}	



	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;
		delete [] elem_error;

	delete grad_1;
	delete grad_2;
	delete dp_1;
	delete dp_2;
	delete mass_matrix;
	delete rhs_1;
	delete rhs_2;
	delete solver_1;
	delete solver_2;

return refine;


}







