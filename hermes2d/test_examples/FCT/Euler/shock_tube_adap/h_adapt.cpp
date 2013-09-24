int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}


	template<typename Scalar>
void calc_elem_error(SpaceSharedPtr<Scalar> space, MeshFunctionSharedPtr<Scalar> sln,MeshFunctionSharedPtr<Scalar> R_h_1,MeshFunctionSharedPtr<Scalar> R_h_2,  HPAdapt* adapt,double h_min, double h_max, 	int* elements_to_refine,	int* no_of_refinement_steps,double* elem_error)
{		

	int ndof = space->get_num_dofs();
	Element* e = NULL;



//-------------- aus reg_estimator------------------------------
	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(sln);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(sln);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);
	CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ; 
	SimpleVector<double> * rhs_1 = new SimpleVector<double>(ndof);
	SimpleVector<double> * rhs_2 = new SimpleVector<double>(ndof);
	dp_1->assemble(mass_matrix,rhs_1); 
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

		int old_steps;

		for_all_active_elements(e, space->get_mesh()){				
			int i = 1;
			old_steps = no_of_refinement_steps[e->id];
			if(elem_error[e->id] >tol_h){ 
										elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; // i++;
					while((elem_error[e->id]> (tol_h+i*std_dev_h))&&(i<5)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			/*}else	if((elem_error[e->id] <EPS)&&(elements_to_refine[e->id] >1)){ 
										elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]< (EPS/i*100))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}*/
			}else {
				if(elements_to_refine[e->id] != 1)	elements_to_refine[e->id] = 0;
			}

			//if(no_of_refinement_steps[e->id]>2) no_of_refinement_steps[e->id]=2;			
			if((elements_to_refine[e->id] == 1)&&(old_steps>no_of_refinement_steps[e->id])) no_of_refinement_steps[e->id]= old_steps;			
			if((elements_to_refine[e->id] == 4)&&(old_steps<no_of_refinement_steps[e->id])&&(old_steps>0)) no_of_refinement_steps[e->id]= old_steps;
			
		}	
		
		
	int max_id =space->get_mesh()->get_max_element_id();
		bool visited[max_id];
		for(int i=0;i<max_id; i++) visited[i] =false;
		int start_id,steps, e_1, e_2,id;
		Element* elem_start_1, *elem_start_2;		

	//Refine all Elements in y-direction
		for_all_active_elements(e, space->get_mesh()) 
		{
			if((elements_to_refine[e->id]==1)&&(visited[e->id]==false))
			{ 
					elem_start_1 = e; start_id = e->id; 
					elem_start_2 = e;
					steps = no_of_refinement_steps[elem_start_1->id];

					if(elem_start_1->vn[0]->y == elem_start_1->vn[1]->y){ 
								e_1 = 0;
								if(elem_start_1->vn[2]->y == elem_start_1->vn[3]->y)
									 e_2 =2;
								else Hermes::Exceptions::Exception("calc_elem_error:falsche y-Koord");
					}else if(elem_start_1->vn[1]->y == elem_start_1->vn[2]->y){
									 e_1 =1;
									if(elem_start_1->vn[3]->y == elem_start_1->vn[0]->y)
										 e_2 =3;
									else Hermes::Exceptions::Exception("calc_elem_error:falsche y-Koord");						 
					}else Hermes::Exceptions::Exception("calc_elem_error:keine zweit knoten selbe y-Koord.");
			
					if( elem_start_1->get_neighbor(e_1)!=NULL) 
					if(no_of_refinement_steps[elem_start_1->get_neighbor(e_1)->id]>steps){ 
								steps = no_of_refinement_steps[elem_start_1->get_neighbor(e_1)->id];
								no_of_refinement_steps[start_id]= steps;
					}
					if( elem_start_2->get_neighbor(e_2)!=NULL) 
					if(no_of_refinement_steps[elem_start_2->get_neighbor(e_2)->id]>steps){ 
							steps = no_of_refinement_steps[elem_start_2->get_neighbor(e_2)->id];
							no_of_refinement_steps[start_id]= steps;
							}

					while((elem_start_1 = elem_start_1->get_neighbor(e_1))!=NULL){
								id = elem_start_1->id;	
								if(visited[id]==false){
									elements_to_refine[id]=1;
									no_of_refinement_steps[id] = steps;
									visited[id] =true;				
								}				
						}
					while((elem_start_2 = elem_start_2->get_neighbor(e_2))!=NULL){ 
								id = elem_start_2->id;	
								if(visited[id]==false){
									elements_to_refine[id]=1;
									no_of_refinement_steps[id] = steps;
									visited[id] =true;				
								}				
					}
			}
		}
	
	//coarse all elements in y-direction		
/*
		for_all_active_elements(e, space->get_mesh()) 
		{
			if((elements_to_refine[e->id]==4)&&(visited[e->id]==false))
			{ 
					elem_start_1 = e; start_id = e->id; 
					elem_start_2 = e;
					steps = no_of_refinement_steps[elem_start_1->id];

					if(elem_start_1->vn[0]->y == elem_start_1->vn[1]->y){ 
								e_1 = 0;
								if(elem_start_1->vn[2]->y == elem_start_1->vn[3]->y)
									 e_2 =2;
								else Hermes::Exceptions::Exception("calc_elem_error:falsche y-Koord");
					}else if(elem_start_1->vn[1]->y == elem_start_1->vn[2]->y){
									 e_1 =1;
									if(elem_start_1->vn[3]->y == elem_start_1->vn[0]->y)
										 e_2 =3;
									else Hermes::Exceptions::Exception("calc_elem_error:falsche y-Koord");						 
					}else Hermes::Exceptions::Exception("calc_elem_error:keine zweit knoten selbe y-Koord.");
			
					if( elem_start_1->get_neighbor(e_1)!=NULL) 
					if(no_of_refinement_steps[elem_start_1->get_neighbor(e_1)->id]<steps){ 
								steps = no_of_refinement_steps[elem_start_1->get_neighbor(e_1)->id];
								no_of_refinement_steps[start_id]= steps;
					}
					if( elem_start_2->get_neighbor(e_2)!=NULL) 
					if(no_of_refinement_steps[elem_start_2->get_neighbor(e_2)->id]<steps){ 
							steps = no_of_refinement_steps[elem_start_2->get_neighbor(e_2)->id];
							no_of_refinement_steps[start_id]= steps;
							}

					while((elem_start_1 = elem_start_1->get_neighbor(e_1))!=NULL){
								id = elem_start_1->id;	
								if(visited[id]==false){
									elements_to_refine[id]=4;
									no_of_refinement_steps[id] = steps;
									visited[id] =true;				
								}				
						}
					while((elem_start_2 = elem_start_2->get_neighbor(e_2))!=NULL){ 
								id = elem_start_2->id;	
								if(visited[id]==false){
									elements_to_refine[id]=4;
									no_of_refinement_steps[id] = steps;
									visited[id] =true;				
								}				
					}
			}
		}
	*/	
		

	delete grad_1;
	delete grad_2;
	delete dp_1;
	delete dp_2;
	delete mass_matrix;
	delete rhs_1;
	delete rhs_2;
	delete solver_1;
	delete solver_2;



}







