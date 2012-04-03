int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}


template<typename Scalar>
 Space<Scalar>* construct_p_refined_space(Space<Scalar>* coarse, int* elements_to_refine)
{
  _F_
	bool changed = false;
	int order, v_ord, h_ord;
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  Space<Scalar>* space = coarse->dup(ref_mesh, 0);
		Element* e = NULL;
	for_all_active_elements(e, space->get_mesh()){
		if(elements_to_refine[e->id]!=1){
			if(e->is_triangle()==true){
				order = space->get_element_order(e->id); 
				if(order >=1){
					if(order<P_MAX) space->set_element_order_internal(e->id, order+1);
					changed = true;
				}
			}else if(e->is_quad()==true){
				v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
				h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
				if((h_ord<P_MAX)&&(v_ord<P_MAX)) order = H2D_MAKE_QUAD_ORDER(h_ord+1, v_ord+1);
				else order = H2D_MAKE_QUAD_ORDER(P_MAX, P_MAX); 
					 space->set_element_order_internal(e->id, order);
					changed = true;
				}		
		}
	}
	space->assign_dofs();
  return space;
}

template<typename Scalar>
 Space<Scalar>* construct_p1_space(Space<Scalar>* coarse)
{
  _F_
	bool changed = false;
	int order, v_ord, h_ord;
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  Space<Scalar>* space = coarse->dup(ref_mesh, 0);
/*		Element* e = NULL;
	for_all_active_elements(e, space->get_mesh()){
			if(e->is_triangle()==true){
				order = space->get_element_order(e->id); 
				if(order >1){
					 space->set_element_order_internal(e->id, 1);
					changed = true;
				}
			}else if(e->is_quad()==true){
				v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
				h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
				if((h_ord>1)||(v_ord>1)){
					 order = H2D_MAKE_QUAD_ORDER(1, 1);			
					 space->set_element_order_internal(e->id, order);
					changed = true;
				}
			}			
		}space->assign_dofs();*/
	space->set_uniform_order(1);	
  return space;
}

template<typename Scalar>
 Space<Scalar>* construct_p_refined_space_with_coarse_mesh(Space<Scalar>* coarse, int* elements_to_refine)
{
  _F_
	bool changed = false;
	int order, v_ord, h_ord;
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  Space<Scalar>* space = coarse->dup(ref_mesh, 0);
		Element* e = NULL;
	for_all_active_elements(e, space->get_mesh()){
		if(elements_to_refine[e->id] !=1) {			//p erhoehen  ->kein FCT
			if(e->is_triangle()==true){
				order = space->get_element_order(e->id); 
				if(order >=1){
					if(order<P_MAX) space->set_element_order_internal(e->id, order+1);
					changed = true;
				}
			}else if(e->is_quad()==true){
				v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
				h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
				if((h_ord<P_MAX)&&(v_ord<P_MAX)) order = H2D_MAKE_QUAD_ORDER(h_ord+1, v_ord+1);
				else order = H2D_MAKE_QUAD_ORDER(P_MAX, P_MAX); 
					 space->set_element_order_internal(e->id, order);
					changed = true;
				}
			}	
		}
	//space->unrefine_all_mesh_elements();
	space->assign_dofs();
  return space;
}
template<typename Scalar>
 Space<Scalar>* construct_space_with_coarse_mesh(Space<Scalar>* coarse, int* elements_to_refine)
{
  _F_
	bool changed = false;
	int order, v_ord, h_ord;
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  Space<Scalar>* space = coarse->dup(ref_mesh, 0);
	space->unrefine_all_mesh_elements();
	space->assign_dofs();
  return space;
}

template<typename Scalar>
bool init_adapt(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* ref_sln, Solution<Scalar>* R_h_1, Solution<Scalar>* R_h_2,CustomWeakFormMassmatrix* massmatrix,HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,  std::list<int>* list,std::list<int>* neighbor, double h_min, double h_max)
{
	if(!list->empty()) list->clear();
	if(!neighbor->empty()) neighbor->clear();
	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	
	int id;
	std::list<int>::iterator it;
	Node* vn;	
	bool in_list = false;
	Element* ne = NULL;
	bool hanging_node =false;

	int* elements_to_refine = new int[space->get_mesh()->get_max_element_id()];   
				int* no_of_refinement_steps = new int[space->get_mesh()->get_max_element_id()];		
	for(int i = 0; i < space->get_mesh()->get_max_element_id(); i++){ 
										elements_to_refine[i]= 0;										
										no_of_refinement_steps[i]=0;
		}
		int dof = space->get_num_dofs();
		int* smooth_fct_in_elem = new int[space->get_mesh()->get_max_element_id()];
		int* smooth_dx_in_elem = new int[space->get_mesh()->get_max_element_id()];
		int* smooth_dy_in_elem = new int[space->get_mesh()->get_max_element_id()];
		int* smooth_elem = new int[space->get_mesh()->get_max_element_id()];
		int* smooth_dof = new int[space->get_num_dofs()];

		Scalar* coeff = new Scalar[dof];Scalar* coeff_2 = new Scalar[dof];
		Scalar* P_plus = new Scalar[dof]; Scalar* P_minus = new Scalar[dof];
		Scalar* Q_plus = new Scalar[dof]; Scalar* Q_minus = new Scalar[dof];	
		Scalar* R_plus = new Scalar[dof]; Scalar* R_minus = new Scalar[dof];	
		DiscreteProblem<Scalar>* dp_mass = new DiscreteProblem<Scalar>(massmatrix, space);
		UMFPackMatrix<Scalar>* mass_matrix = new UMFPackMatrix<Scalar>; 
		dp_mass->assemble(mass_matrix); 	
		UMFPackMatrix<Scalar>* lumped_matrix = massLumping(dof_list,dof_list_2,mass_matrix);
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
		Lumped_Projection::project_lumped(space, sln, coeff, matrix_solver, lumped_matrix);
		Solution<Scalar>::vector_to_solution(coeff, space, ref_sln);
	smoothness_indicator(space,ref_sln,R_h_1,R_h_2, smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem,smooth_elem,smooth_dof,al);
		OGProjection<Scalar>::project_global(space,sln, coeff_2, matrix_solver, HERMES_L2_NORM);
		lumped_flux_limiter(dof_list,dof_list_2,mass_matrix, lumped_matrix, coeff, coeff_2,
								P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus );
		Solution<Scalar>::vector_to_solution(coeff, space, ref_sln);


	Adapt<Scalar>* adaptivity = new Adapt<Scalar>(space, HERMES_H1_NORM);

	double err_est = adaptivity->calc_err_est(sln,ref_sln,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);
		double mean_value=0.;
		double std_dev = 0.0;	
			int counter = space->get_mesh()->get_num_active_elements();	
		for_all_active_elements(e, space->get_mesh()){
			mean_value += std::sqrt(adaptivity->get_element_error_squared(0, e->id));
			std_dev += adaptivity->get_element_error_squared(0, e->id);
		}	
			std_dev -= ((mean_value*mean_value)/counter); mean_value/=counter;
		std_dev/=(counter-1);
		std_dev = std::sqrt(std_dev);
			
			double tol = mean_value;
			double tol_2 = mean_value + std_dev;



		for_all_active_elements(e, space->get_mesh()){
				if(std::sqrt(adaptivity->get_element_error_squared(0, e->id))>tol){						
					refine = true; elements_to_refine[e->id] =1;no_of_refinement_steps[e->id]++;
				}else if(std::sqrt(adaptivity->get_element_error_squared(0, e->id))<EPS)	{
					refine = true; elements_to_refine[e->id] =4;no_of_refinement_steps[e->id]++;
				}else if(smooth_elem[e->id]==1)		{
								refine = true; elements_to_refine[e->id] =4;no_of_refinement_steps[e->id]++;
				}
		}

	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

			delete [] P_plus; 
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete [] coeff;
			delete [] coeff_2;
			delete dp_mass;
			delete lumped_matrix; 
			delete mass_matrix; 
			delete adaptivity;
			delete [] elements_to_refine;
			delete [] no_of_refinement_steps;

		delete [] smooth_fct_in_elem;
		delete [] smooth_dx_in_elem;
		delete [] smooth_dy_in_elem;
		delete [] smooth_elem;
		delete [] smooth_dof;

	return refine;

}

	template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* p1_sln,Solution<Scalar>* u_2h,Solution<Scalar>* u_2h_2p,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,  ResidualForm* residual,CustomWeakFormMassmatrix* massmatrix, HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,  std::list<int>* list,std::list<int>* neighbor, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps)
{		

	if(!list->empty()) list->clear();
	if(!neighbor->empty()) neighbor->clear();

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
		double tol_h = mean_value_h+std_dev_h;

if(ps==1){
		for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if(elem_error[e->id] >tol_h){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; i++;
					while((elem_error[e->id]> (tol_h+i*std_dev_h))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}
		}
}else{	
	for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if(elem_error[e->id] <EPS){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; i++;
					while((elem_error[e->id]< (EPS/i*1000))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}
		}
	


}
		delete [] elem_error;

	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

return refine;


}











/*

template<typename Scalar>
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* p1_sln,Solution<Scalar>* u_2h,Solution<Scalar>* u_2h_2p,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2,  ResidualForm* residual,CustomWeakFormMassmatrix* massmatrix, HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,  std::list<int>* list,std::list<int>* neighbor, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps)
{		

	if(!list->empty()) list->clear();
	if(!neighbor->empty()) neighbor->clear();

	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	

	int id;
	std::list<int>::iterator it;
	Node* vn;	
	bool in_list = false;
	Element* ne = NULL;
	bool hanging_node =false;



//-----------h-Adapt

	if(ps==1){

			//----------h-Indicator for Refinement
		Scalar* u_h =NULL;
		Space<Scalar>* ref_space_h = construct_p1_space(space);
		DiscreteProblem<Scalar>* dp_h = new DiscreteProblem<Scalar>(residual, ref_space_h);
		int dof= ref_space_h->get_num_dofs();
		UMFPackMatrix<Scalar>* matrix_h = new UMFPackMatrix<Scalar>;
		UMFPackVector<Scalar>* rhs_h = new UMFPackVector<Scalar>(dof); 
		dp_h->assemble(matrix_h, rhs_h);
				  // Solve the linear system and if successful, obtain the solution. 
		UMFPackLinearSolver<Scalar>* solver_h = new UMFPackLinearSolver<Scalar>(matrix_h,rhs_h);	
		if(solver_h->solve()){ 
			u_h = solver_h->get_sln_vector();  
			Solution<Scalar>::vector_to_solution(u_h, ref_space_h, u_2h);	
		}else error ("Matrix solver failed.\n");
		Adapt<Scalar>* adaptivity_e_h = new Adapt<Scalar>(ref_space_h, HERMES_L2_NORM);
		
		Scalar* zeros = new Scalar[dof]; for(int i =0; i<dof; i++) zeros[i] = 0.;
		Solution<Scalar>* zero_sln = new Solution<Scalar>(ref_space_h, zeros);

double err_est_e_h = adaptivity_e_h->calc_err_est(zero_sln,u_2h,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);

		Scalar* zeros_2 = new Scalar[space->get_num_dofs()]; for(int i =0; i<space->get_num_dofs(); i++) zeros_2[i] = 0.;
		Solution<Scalar>* zero_sln_2 = new Solution<Scalar>(space, zeros_2);
			Adapt<Scalar>* adaptivity_h = new Adapt<Scalar>(space ,HERMES_H1_SEMINORM);
		double err_est_h = adaptivity_h->calc_err_est(zero_sln_2,sln,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);

		double mean_value_e_h=0.;
		double std_dev_e_h = 0.0;
		double mean_value_h=0.;
		double std_dev_h = 0.0;
		int counter = space->get_mesh()->get_num_active_elements();	
		for_all_active_elements(e, space->get_mesh()){
			if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){
				mean_value_e_h += std::sqrt(adaptivity_e_h->get_element_error_squared(0, e->id));
				std_dev_e_h += adaptivity_e_h->get_element_error_squared(0, e->id);		
				mean_value_h += std::sqrt(adaptivity_h->get_element_error_squared(0, e->id));
				std_dev_h += adaptivity_h->get_element_error_squared(0, e->id);	
			}
		}			

	std_dev_h -= ((mean_value_h*mean_value_h)/counter); mean_value_h/=counter;
		std_dev_h/=(counter-1);
		std_dev_h = std::sqrt(std_dev_h);
	double tol_h = mean_value_h;
	double tol_h_min = mean_value_h - 0.5*std_dev_h;


		std_dev_e_h -= ((mean_value_e_h*mean_value_e_h)/counter); mean_value_e_h/=counter;
		std_dev_e_h/=(counter-1);
		std_dev_e_h = std::sqrt(std_dev_e_h);			


		double tol_e_h = mean_value_e_h +std_dev_e_h;

		double tol_e_h_2 = mean_value_e_h - std_dev_e_h;
		if(tol_e_h_2<0.) tol_e_h_2 =EPS;

	
		for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if((std::sqrt(adaptivity_e_h->get_element_error_squared(0, e->id)) >tol_e_h)){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((std::sqrt(adaptivity_e_h->get_element_error_squared(0, e->id))> (tol_e_h+i*std_dev_e_h))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else if((std::sqrt(adaptivity_h->get_element_error_squared(0, e->id)) >tol_h)){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((std::sqrt(adaptivity_h->get_element_error_squared(0, e->id))> (tol_e_h+i*std_dev_h))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}
		}

		delete matrix_h;
		delete rhs_h;
		delete solver_h;
		delete dp_h;
		delete ref_space_h->get_mesh();
		delete ref_space_h;
		delete adaptivity_e_h;
			delete zero_sln;
	delete [] zeros;
		delete adaptivity_h;
	delete zero_sln_2;
	delete [] zeros_2;
		



	}else{
			//----------h-Indicator for Coarsening

			Space<Scalar>* ref_space_h = construct_space_with_coarse_mesh(space, elements_to_refine);
			Adapt<Scalar>* adaptivity_h = new Adapt<Scalar>(space, HERMES_L2_NORM);
				int dof_h = ref_space_h->get_num_dofs();
			Scalar* coeff = new Scalar[dof_h];

			AsmList<Scalar> dof_list_h;
			AsmList<Scalar> dof_list_h_2;
			p1_list_fast(ref_space_h, &dof_list_h,&dof_list_h_2, al);
			Scalar* coeff_3 = new Scalar[dof_h];
			Scalar* P_plus = new Scalar[dof_h]; Scalar* P_minus = new Scalar[dof_h];
			Scalar* Q_plus = new Scalar[dof_h]; Scalar* Q_minus = new Scalar[dof_h];	
			Scalar* R_plus = new Scalar[dof_h]; Scalar* R_minus = new Scalar[dof_h];	
			DiscreteProblem<Scalar>* dp_mass = new DiscreteProblem<Scalar>(massmatrix, ref_space_h);
			UMFPackMatrix<Scalar>* mass_matrix = new UMFPackMatrix<Scalar>; 
			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<Scalar>* lumped_matrix = massLumping(&dof_list_h,&dof_list_h_2,mass_matrix);
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
			Lumped_Projection::project_lumped(ref_space_h, sln, coeff, matrix_solver, lumped_matrix);
			OGProjection<Scalar>::project_global(ref_space_h,sln, coeff_3, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(&dof_list_h,&dof_list_h_2,mass_matrix, lumped_matrix, coeff, coeff_3,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus );
			Solution<Scalar>::vector_to_solution(coeff, ref_space_h, u_2h);


			//OGProjection::project_global(ref_space_h,sln, coeff, matrix_solver, HERMES_L2_NORM);
			//Solution::vector_to_solution(coeff, ref_space_h, u_2h);
			double	err_est = adaptivity_h->calc_err_est(sln,u_2h,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS) ;

			double mean_value_h=0.;
			double std_dev_h = 0.0;	
			int counter = space->get_mesh()->get_num_active_elements();	
			for_all_active_elements(e, space->get_mesh()){
				mean_value_h += std::sqrt(adaptivity_h->get_element_error_squared(0, e->id));
				std_dev_h += adaptivity_h->get_element_error_squared(0, e->id);
			}			
			std_dev_h -= ((mean_value_h*mean_value_h)/counter); mean_value_h/=counter;
			std_dev_h/=(counter-1);
			std_dev_h = std::sqrt(std_dev_h);
			
			double tol_h = mean_value_h-std_dev_h;			
			if(tol_h<0.) tol_h = EPS;	
			//double tol_h =EPS;			

			for_all_active_elements(e, space->get_mesh()){	
				int i =0;
				if((std::sqrt(adaptivity_h->get_element_error_squared(0, e->id))< tol_h)){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; i++;
						while((std::sqrt(adaptivity_h->get_element_error_squared(0, e->id))<(tol_h/(i*1000.)))&&(i<2)){
							 													no_of_refinement_steps[e->id]++; i++;
							}
				}
			}

			delete [] P_plus; 
			delete [] P_minus;
			delete [] Q_plus;
			delete [] Q_minus;
			delete [] R_plus;
			delete [] R_minus;
			delete [] coeff_3;
			delete dp_mass;
			delete lumped_matrix; 
			delete mass_matrix; 



			delete ref_space_h->get_mesh();
			delete ref_space_h;
			delete adaptivity_h;
			delete [] coeff;

	}




	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);

/*
//------------------Elemente fuer FCT speichern
	 bool p2_neighbor =false;
	double elem_diag =0; Element* elem_neigh;
	for_all_active_elements(e, space->get_mesh()){  
		elem_diag=e->get_diameter();
		if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
							 	elem_neigh = e->get_neighbor(iv);
								if(elem_neigh!=NULL){ 
										id = elem_neigh->id;	
										if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
											p2_neighbor =true; neighbor->push_back(e->id);
											break;
										}else if(elem_diag != elem_neigh->get_diameter()){
													// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neighbor =true; neighbor->push_back(e->id);
											break;
										}
								}
						}
						if(p2_neighbor==false){list->push_back(e->id);					
						}  // Elemente fuer FCT
						else {p2_neighbor =false;				
						}
		}else if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(e->id)==2)){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
							 	elem_neigh = e->get_neighbor(iv);
								if(elem_neigh!=NULL){ 
										id = elem_neigh->id;	
										if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(id)==1)){
											p2_neighbor =true; neighbor->push_back(e->id);
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

}*/

