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
	space->unrefine_all_mesh_elements();
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
/*bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* u_prev_time, Solution<Scalar>* sln,Solution<Scalar>* u_2h,Solution<Scalar>* u_2h_2p,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2, CustomWeakFormMassmatrix* massmatrix, HPAdapt* adapt, AsmList<Scalar>* dof_list,AsmList<Scalar>* dof_list_2, AsmList<Scalar>* al,  std::list<int>* list,std::list<int>* neighbor, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps, int* smooth_elem)*/
bool h_p_adap(Space<Scalar>* space,Solution<Scalar>* sln,Solution<Scalar>* R_h_1,Solution<Scalar>* R_h_2, HPAdapt* adapt, std::list<int>* list,std::list<int>* neighbor, int* elements_to_refine,int* no_of_refinement_steps,double h_min, double h_max, int ts, int ps, int* smooth_elem)
{		

	if(!list->empty()) list->clear();
	if(!neighbor->empty()) neighbor->clear();

	Element* e = NULL;
	bool refine = false;  //->verfeinern oder vergroebern?	
int id;



//-----------h-Adapt	

		double* elem_error = new double[space->get_mesh()->get_max_element_id()];
		calc_z_z_error( space, sln, R_h_1, R_h_2, elem_error);
		double mean_value_z=0.;
		double std_dev_z = 0.0;
		int counter = space->get_mesh()->get_num_active_elements();
		for_all_active_elements(e, space->get_mesh()){
				mean_value_z += elem_error[e->id];
				std_dev_z += elem_error[e->id]*elem_error[e->id];			
		}	

		std_dev_z -= ((mean_value_z*mean_value_z)/counter); mean_value_z/=counter;
		std_dev_z/=(counter-1);
		std_dev_z = std::sqrt(std_dev_z);
		double tol_z = mean_value_z+std_dev_z;

		for_all_active_elements(e, space->get_mesh()){	
			int i = 1;	
			if(elem_error[e->id] >tol_z){ 
										refine = true; elements_to_refine[e->id] = 1; no_of_refinement_steps[e->id]++; 
					while((elem_error[e->id]> (tol_z+i*std_dev_z))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else	if(smooth_elem[e->id]==1){
							refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;
			}else if(elem_error[e->id] <EPS){ 
										refine = true; elements_to_refine[e->id] = 4; no_of_refinement_steps[e->id]++; i++;
					while((elem_error[e->id]< (EPS/i*1000))&&(i<2)){
					 													no_of_refinement_steps[e->id]++; i++;
					}
			}else{ elements_to_refine[e->id] =1;no_of_refinement_steps[e->id]++;}
		}

		//-----------p-Adapt

		/*	Space<Scalar>* ref_space_p = construct_p_refined_space_with_coarse_mesh(space, elements_to_refine);
			int dof_p = ref_space_p->get_num_dofs();
			Scalar* coeff_2 = new Scalar[dof_p];
			//OGProjection::project_global(ref_space_p,sln, coeff_2, matrix_solver, HERMES_L2_NORM);
			//Solution::vector_to_solution(coeff_2, ref_space_p, u_2h_2p);

			AsmList<Scalar> dof_list_p;
			AsmList<Scalar> dof_list_p_2;
			p1_list_fast(ref_space_p, &dof_list_p,&dof_list_p_2, al);
			Scalar* coeff_3 = new Scalar[dof_p];
			Scalar* P_plus = new Scalar[dof_p]; Scalar* P_minus = new Scalar[dof_p];
			Scalar* Q_plus = new Scalar[dof_p]; Scalar* Q_minus = new Scalar[dof_p];	
			Scalar* R_plus = new Scalar[dof_p]; Scalar* R_minus = new Scalar[dof_p];	
			DiscreteProblem<Scalar>* dp_mass = new DiscreteProblem<Scalar>(massmatrix, ref_space_p);
			UMFPackMatrix<Scalar>* mass_matrix = new UMFPackMatrix<Scalar>; 
			dp_mass->assemble(mass_matrix); 	
			UMFPackMatrix<Scalar>* lumped_matrix = massLumping(&dof_list_p,&dof_list_p_2,mass_matrix);
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
			Lumped_Projection::project_lumped(ref_space_p, sln, coeff_2, matrix_solver, lumped_matrix);
			OGProjection<Scalar>::project_global(ref_space_p,sln, coeff_3, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(&dof_list_p,&dof_list_p_2,mass_matrix, lumped_matrix, coeff_2, coeff_3,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus );
			Solution<Scalar>::vector_to_solution(coeff_2, ref_space_p, u_2h_2p);

			Adapt<Scalar>* adaptivity_p = new Adapt<Scalar>(space, HERMES_H1_NORM);	
			double err_est = adaptivity_p->calc_err_est(sln,u_2h_2p,true,HERMES_TOTAL_ERROR_ABS|HERMES_ELEMENT_ERROR_ABS);


			double mean_value=0.;
			double std_dev = 0.0;	
			int counter_2 = 0;		
			for_all_active_elements(e, space->get_mesh()){
				if(elements_to_refine[e->id] != 1){
					mean_value += std::sqrt(adaptivity_p->get_element_error_squared(0, e->id));
					std_dev += adaptivity_p->get_element_error_squared(0, e->id);
					counter_2++;
				}		
			}	
			std_dev -= ((mean_value*mean_value)/counter_2); mean_value/=counter_2;
			std_dev/=(counter_2-1);
			std_dev = std::sqrt(std_dev);
					
		double tol = mean_value;
		double tol_min = mean_value - 0.25*std_dev;


		for_all_active_elements(e, space->get_mesh()){	
		//	if(elements_to_refine[e->id] != 1){
			//if((elements_to_refine[e->id] != 1)	&&(std::sqrt(adaptivity_p->get_element_error_squared(0, e->id))>tol))	
					if(smooth_elem[e->id]==1){
							refine = true; elements_to_refine[e->id] =2; no_of_refinement_steps[e->id]++;
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

			delete adaptivity_p;
			delete [] coeff_2;
			delete ref_space_p->get_mesh();
			delete ref_space_p;
*/

		

delete [] elem_error;


	Element* elem_neigh=NULL;
	bool p2_neighbor = false;
	for_all_active_elements(e, space->get_mesh()){
			if(elements_to_refine[e->id]==2){//elements_to_refine[e->id]=1;
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
					if(p2_neighbor==false) elements_to_refine[e->id]=0;
					else p2_neighbor = false;			
				}
	}

	

	/*for_all_active_elements(e, space->get_mesh()){
			if(elements_to_refine[e->id]==2){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 
								id = elem_neigh->id;	
								if(elements_to_refine[id]==0){elements_to_refine[e->id]=1;no_of_refinement_steps[e->id]=1;
								}
							}
					}
		
				}
	}*/




	if(refine==true) refine = adapt->adapt(elements_to_refine,no_of_refinement_steps,P_MAX, h_min,h_max,NDOF_STOP);




//------------------Elemente fuer FCT speichern
	 p2_neighbor =false;
	double elem_diag =0; 
	for_all_active_elements(e, space->get_mesh()){  
		elem_diag=e->get_diameter();
		if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
							 	elem_neigh = e->get_neighbor(iv);
								if(elem_neigh!=NULL){ 
										id = elem_neigh->id;	
										if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2)){
											p2_neighbor =true; //neighbor->push_back(e->id);
											break;
										}else if(elem_diag != elem_neigh->get_diameter()){
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
		}
	}

	return refine;

}

