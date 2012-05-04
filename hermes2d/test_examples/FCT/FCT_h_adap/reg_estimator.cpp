
//Fuer approx. Loesungen: Solution<double>* sln

double linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

			//Mittelpunkt des Referenzelements (Einheitsquadrat)	
			double x_c_ref = 0.;
			double y_c_ref = 0.; 

			double u_h_x_c = sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
			double u_h_hat = u_h_x_c + R_h_1->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(x_i-x_c)+R_h_2->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(y_i-y_c);
	
	return u_h_hat;
}

double linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 

		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);

		double u_h_hat = d_u_h_x_c + R_h_1->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)
															 + R_h_2->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(y_i-y_c);
return u_h_hat;

}

double linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 

		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);

		double u_h_hat = d_u_h_x_c + R_h_1->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(x_i-x_c)
																+R_h_2->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);
return u_h_hat;

}

void smoothness_indicator(Space<double>* space,Solution<double>* sln,Solution<double>* R_h_1,Solution<double>* R_h_2, int* smooth_elem_patch, int* smooth_dof, AsmList<double>* al, UMFPackMatrix<double> * mass_matrix ){

	if(sln==NULL) error("smoothness_indicator: sln=NULL");
	if(space==NULL) error("smoothness_indicator: space=NULL");

	int ndof = space->get_num_dofs();


	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(sln);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(sln);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);
	UMFPackVector<double> * rhs_1 = new UMFPackVector<double>(ndof);
	UMFPackVector<double> * rhs_2 = new UMFPackVector<double>(ndof);
	dp_1->assemble(rhs_1); 
	dp_2->assemble(rhs_2);
	UMFPackLinearSolver<double> * solver_1 = new UMFPackLinearSolver<double> (mass_matrix,rhs_1);
	if(solver_1->solve()){ 				
		Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , space, R_h_1);	
	}else error ("Matrix solver failed.\n");
	UMFPackLinearSolver<double> * solver_2 = new UMFPackLinearSolver<double> (mass_matrix,rhs_2);	
	if(solver_2->solve()){ 				
		Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , space, R_h_2);	
	}else error ("Matrix solver failed.\n");



//---------------Elemente mit Vertex(dof) bestimmen
	std::list<int> dof_elem_list[ndof];	
		double* u_c = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dx = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dy = new double[space->get_mesh()->get_max_element_id()];
				//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
		Element* e =NULL; 
		int* index = new int[4];
		for(int i =0;i<4;i++)  index[i] =  space->get_shapeset()->get_vertex_index(i,HERMES_MODE_QUAD);

	for_all_active_elements(e, space->get_mesh()){
			u_c[e->id]= sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
			d_u_c_dx[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);
			d_u_c_dy[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);
			space->get_element_assembly_list(e, al);
			for (unsigned int iv = 0; iv < e->get_nvert(); iv++){ 
							for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
    							if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0)){
											dof_elem_list[al->get_dof()[j]].push_back(e->id); 	
									}
							}
			}	

	}

//------------------------------------------------

	Node* vn=NULL;
	double x_c =0.; double y_c =0.; double u_h_x_c =0.;
	double* x = new double[4]; double* y = new double[4];
	double u_i; double u_dx ,u_dy;
	double u_min, u_max, u_min_dx, u_min_dy, u_max_dx, u_max_dy;
	std::list<int>::iterator elem_id; 
	int smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem;
	double epsilon = EPS; 

	bool non_smooth = false; bool non_smooth_dx = false;bool non_smooth_dy = false;

	for_all_active_elements(e, space->get_mesh()){
			non_smooth = false; non_smooth_dx = false; non_smooth_dy = false;
			smooth_fct_in_elem =1; //erstmal als smooth annehmen	
			smooth_dx_in_elem =1; 
			smooth_dy_in_elem =1; 
			smooth_elem_patch[e->id]=0;
			space->get_element_assembly_list(e, al);		
				x_c = 0.; y_c = 0.;
			for (unsigned int iv = 0; iv < e->get_nvert(); iv++){ 		 
				 vn = e->vn[iv];
					x[iv]=vn->x; y[iv]=vn->y;  //x/y-Koordinaten der Knotenpunkte
					x_c+= (vn->x); y_c+= (vn->y);   // Mittelpunkte bestimmen
			}
			x_c/=e->get_nvert(); y_c/=e->get_nvert();
			for (unsigned int iv = 0; iv < e->get_nvert(); iv++){ 			
				for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
						if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0)){
								u_i = linear_approx(e, x[iv], y[iv],x_c, y_c,sln, R_h_1,R_h_2);	
								u_dx =	linear_approx_dx(e, x[iv], y[iv],x_c, y_c,sln, R_h_1,R_h_2);
								u_dy =	linear_approx_dy(e, x[iv], y[iv],x_c, y_c,sln, R_h_1,R_h_2);	
								u_min = u_c[e->id]; u_max =u_c[e->id];	
								u_min_dx= d_u_c_dx[e->id];	u_max_dx= d_u_c_dx[e->id];	
								u_min_dy= d_u_c_dy[e->id];	u_max_dy= d_u_c_dy[e->id];								
								for(elem_id=dof_elem_list[al->get_dof()[j]].begin();elem_id!=dof_elem_list[al->get_dof()[j]].end();elem_id++){	
									if(u_c[*elem_id]>u_max) u_max =	u_c[*elem_id];
									else if(u_c[*elem_id]<u_min) u_min =	u_c[*elem_id];	
									if(d_u_c_dx[*elem_id]>u_max_dx) u_max_dx =	d_u_c_dx[*elem_id];
									else if(d_u_c_dx[*elem_id]<u_min_dx) u_min_dx =	d_u_c_dx[*elem_id];	
									if(d_u_c_dy[*elem_id]>u_max_dy) u_max_dy =	d_u_c_dy[*elem_id];
									else if(d_u_c_dy[*elem_id]<u_min_dy) u_min_dy =	d_u_c_dy[*elem_id];								
								}
								if((u_i>u_min+epsilon)&&(u_i<u_max-epsilon)) non_smooth = false;
								else non_smooth = true;
								if((u_dx>u_min_dx+epsilon)&&(u_dx<u_max_dx-epsilon)) non_smooth_dx = false;
								else non_smooth_dx = true;
								if((u_dy>u_min_dy+epsilon)&&(u_dy<u_max_dy-epsilon)) non_smooth_dy = false;
								else non_smooth_dy = true;			
								break;
						}
				}
				if(non_smooth == true) smooth_fct_in_elem=0; 
				if(non_smooth_dx == true) smooth_dx_in_elem=0; 
				if(non_smooth_dy == true) smooth_dy_in_elem=0; 
				if((non_smooth == true)&&(non_smooth_dx == true)&&(non_smooth_dy == true)) break;
			}
			if(max(smooth_fct_in_elem, min(smooth_dx_in_elem,smooth_dy_in_elem))==1){
						smooth_elem_patch[e->id]=1;
			}
	}


	
for(int i =0; i<ndof;i++){
			non_smooth = false; 
			smooth_dof[i]=0;
			for(elem_id=dof_elem_list[i].begin();elem_id!=dof_elem_list[i].end();elem_id++){
					if(smooth_elem_patch[*elem_id]==0){ non_smooth = true; break;}
			}
			if(non_smooth==true) {
					for(elem_id=dof_elem_list[i].begin();elem_id!=dof_elem_list[i].end();elem_id++){
						if(smooth_elem_patch[*elem_id]==1)	smooth_elem_patch[*elem_id]=2;
					}			
			}else{
						smooth_dof[i]=1;
			}
	}

	


  //Clean-up
	delete grad_1;
	delete grad_2;
	delete dp_1;
	delete dp_2;

	delete rhs_1;
	delete rhs_2;
	delete solver_1;
	delete solver_2;
	delete [] u_c;
	delete [] d_u_c_dx ;
	delete [] d_u_c_dy ;
	delete [] index;
	delete [] x ;
	delete [] y;


}


