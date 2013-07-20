int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}


	template<typename Scalar>
bool h_p_adap(SpaceSharedPtr<Scalar> space,	UMFPackMatrix<double> * mass_matrix,MeshFunctionSharedPtr<Scalar> sln,MeshFunctionSharedPtr<Scalar> R_h_1,MeshFunctionSharedPtr<Scalar> R_h_2,  HPAdapt* adapt,double h_min, double h_max,	int* elements_to_refine,	int* no_of_refinement_steps,		double* elem_error)
{		

	int ndof = space->get_num_dofs();
	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	


//-------------- aus reg_estimator------------------------------
	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(sln);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(sln);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);
	//UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ; 
	UMFPackVector<double> * rhs_1 = new UMFPackVector<double>(ndof);
	UMFPackVector<double> * rhs_2 = new UMFPackVector<double>(ndof);
	dp_1->assemble(rhs_1); 
	dp_2->assemble(rhs_2);
	UMFPackLinearMatrixSolver<double> * solver_1 = new UMFPackLinearMatrixSolver<double> (mass_matrix,rhs_1);
  try
  {
   solver_1->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }					
		Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , space, R_h_1);	

	UMFPackLinearMatrixSolver<double> * solver_2 = new UMFPackLinearMatrixSolver<double> (mass_matrix,rhs_2);	
  try
  {
   solver_2->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }				
		Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , space, R_h_2);	

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



	delete grad_1;
	delete grad_2;
	delete dp_1;
	delete dp_2;
	//delete mass_matrix;
	delete rhs_1;
	delete rhs_2;
	delete solver_1;
	delete solver_2;

return refine;


}







