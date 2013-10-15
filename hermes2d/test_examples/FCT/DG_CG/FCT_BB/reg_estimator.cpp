int min(int a, int b){
		if(a>b) return b;
		else a;	
}
int max(int a, int b){
		if(a<b) return b;
		else a;	
}

ElementMode2D mode = HERMES_MODE_TRIANGLE;


//Fuer approx. Loesungen: Solution<double>* sln

double linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

			//Mittelpunkt des Referenzelements (Einheitsquadrat)	
			double x_c_ref = 0.;
			double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}

			double u_h_x_c = sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
			double u_h_hat = u_h_x_c + R_h_1->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(x_i-x_c)+R_h_2->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(y_i-y_c);
	
	return u_h_hat;
}

double linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}

		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);

		double u_h_hat = d_u_h_x_c + R_h_1->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)															 + R_h_1->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);


return u_h_hat;

}

double linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln, Solution<double>* R_h_1,Solution<double>* R_h_2){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}

		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);

		double u_h_hat = d_u_h_x_c + R_h_2->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)															+R_h_2->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);


return u_h_hat;

}



double linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln){

			//Mittelpunkt des Referenzelements (Einheitsquadrat)	
			double x_c_ref = 0.;
			double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}
			double u_h_x_c = sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);

double u_h_hat = u_h_x_c + sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)+sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);
	
	return u_h_hat;
}

double linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}
		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);

		double u_h_hat = d_u_h_x_c + sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 3)*(x_i-x_c)															 + sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 5)*(y_i-y_c);
return u_h_hat;

}

double linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln){

	//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}
		double d_u_h_x_c = sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);

double u_h_hat = d_u_h_x_c + sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 5)*(x_i-x_c)		        															+sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 4)*(y_i-y_c);
return u_h_hat;

}



void smoothness_indicator(SpaceSharedPtr<double> space,MeshFunctionSharedPtr<double> sln_fct,int* smooth_elem_patch, int* smooth_dof, AsmList<double>* al,CSCMatrix<double> * mass_matrix,bool proj){

	if(sln_fct==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: sln=NULL");
	if(space==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: space=NULL");

	int ndof = space->get_num_dofs();

Solution<double>* sln = new Solution<double>; 
sln = static_cast<Solution<double>* > (sln_fct->clone());


//---------------Elemente mit Vertex(dof) bestimmen
	std::list<int> dof_elem_list[ndof];	
		double* u_c = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dx = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dy = new double[space->get_mesh()->get_max_element_id()];
				//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}
		Element* e =NULL; int* index;
	if(mode == HERMES_MODE_QUAD){
		index = new int[4];
		for(int i =0;i<4;i++)  index[i] =  space->get_shapeset()->get_vertex_index(i,mode);
	}else{
		index = new int[3];
		for(int i =0;i<3;i++)  index[i] =  space->get_shapeset()->get_vertex_index(i,mode);
	}
	for_all_active_elements(e, space->get_mesh()){
			u_c[e->id]= sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
			d_u_c_dx[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);
			d_u_c_dy[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);
			space->get_element_assembly_list(e, al);
			//for (unsigned int iv = 0; iv < e->get_nvert(); iv++){ 
							for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
    							//if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0))
									if(al->get_dof()[j]!=-1.0){
											dof_elem_list[al->get_dof()[j]].push_back(e->id); 	
									}
							}
			//}	

	}

//------------------------------------------------

	Node* vn=NULL;
	double x_c =0.; double y_c =0.; double u_h_x_c =0.;
	double* x = new double[4]; double* y = new double[4];
	double u_i; double u_dx ,u_dy;
	double u_min, u_max, u_min_dx, u_min_dy, u_max_dx, u_max_dy;
	std::list<int>::iterator elem_id; 
			 int smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem;
	double epsilon = EPS_smooth; 

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
							u_i = linear_approx(e, x[iv], y[iv],x_c, y_c,sln);	
							u_dx =	linear_approx_dx(e, x[iv], y[iv],x_c, y_c,sln);
							u_dy =	linear_approx_dy(e, x[iv], y[iv],x_c, y_c,sln);	
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
			if(proj==true){
				if(max(smooth_fct_in_elem, min(smooth_dx_in_elem,smooth_dy_in_elem))==1){
						smooth_elem_patch[e->id]=1;  
						e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("s").marker;
				}		
				else{ e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("").marker;}
			}else{
				if(min(smooth_dx_in_elem,smooth_dy_in_elem)==1)	{
						smooth_elem_patch[e->id]=1;  
						e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("s").marker;
				}
				else{ e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("").marker;
				}
			}
	}



for(int i =0; i<ndof;i++){
			non_smooth = false; smooth_dof[i]=0;
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
	delete [] u_c;
	delete [] d_u_c_dx ;
	delete [] d_u_c_dy ;
	delete [] index;
	delete [] x ;
	delete [] y;

delete sln; 
}









void smoothness_indicator(SpaceSharedPtr<double> space,MeshFunctionSharedPtr<double> sln_fct,MeshFunctionSharedPtr<double> R_h_1_fct,MeshFunctionSharedPtr<double> R_h_2_fct, int* smooth_elem_patch, int* smooth_dof, AsmList<double>* al,CSCMatrix<double> * mass_matrix,bool proj){

	if(sln_fct==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: sln=NULL");
	if(space==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: space=NULL");

	int ndof = space->get_num_dofs();


	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(sln_fct);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(sln_fct);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);

	SimpleVector<double> * rhs_1 = new SimpleVector<double>(ndof);
	SimpleVector<double> * rhs_2 = new SimpleVector<double>(ndof);
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
		Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , space, R_h_1_fct);	

	UMFPackLinearMatrixSolver<double> * solver_2 = new UMFPackLinearMatrixSolver<double> (mass_matrix,rhs_2);	
  try
  {
   solver_2->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }					
		Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , space, R_h_2_fct);	


Solution<double>* sln = new Solution<double>; sln = static_cast<Solution<double>* > (sln_fct->clone());
Solution<double>* R_h_1 = new Solution<double>; R_h_1= static_cast<Solution<double>* > (R_h_1_fct->clone());
Solution<double>* R_h_2 = new Solution<double>; R_h_2= static_cast<Solution<double>* > (R_h_2_fct->clone());

//---------------Elemente mit Vertex(dof) bestimmen
	std::list<int> dof_elem_list[ndof];	
		double* u_c = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dx = new double[space->get_mesh()->get_max_element_id()];
		double* d_u_c_dy = new double[space->get_mesh()->get_max_element_id()];
				//Mittelpunkt des Referenzelements (Einheitsquadrat)	
		double x_c_ref = 0.;
		double y_c_ref = 0.; 
	if(mode == HERMES_MODE_TRIANGLE){x_c_ref = -2./3.; y_c_ref = -2./3.;}
		Element* e =NULL; 
		int* index = new int[4];
	if(mode == HERMES_MODE_QUAD){
		index = new int[4];
		for(int i =0;i<4;i++)  index[i] =  space->get_shapeset()->get_vertex_index(i,mode);
	}else{
		index = new int[3];
		for(int i =0;i<3;i++)  index[i] =  space->get_shapeset()->get_vertex_index(i,mode);
	}

	for_all_active_elements(e, space->get_mesh()){
			u_c[e->id]= sln->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
			d_u_c_dx[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);
			d_u_c_dy[e->id]= sln->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);
			space->get_element_assembly_list(e, al);
			//for (unsigned int iv = 0; iv < e->get_nvert(); iv++){ 
							for(unsigned int j = 0; j < al->get_cnt(); j ++){			 
    							//if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0))
									if(al->get_dof()[j]!=-1.0){
											dof_elem_list[al->get_dof()[j]].push_back(e->id); 	
									}
							}
			//}	

	}

//------------------------------------------------

	Node* vn=NULL;
	double x_c =0.; double y_c =0.; double u_h_x_c =0.;
	double* x = new double[4]; double* y = new double[4];
	double u_i; double u_dx ,u_dy;
	double u_min, u_max, u_min_dx, u_min_dy, u_max_dx, u_max_dy;
	std::list<int>::iterator elem_id; 
			 int smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem;
	double epsilon = EPS_smooth; 

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
			if(proj==true){
				if(max(smooth_fct_in_elem, min(smooth_dx_in_elem,smooth_dy_in_elem))==1){
						smooth_elem_patch[e->id]=1;  
						e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("s").marker;
				}		
				else{ e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("").marker;}
			}else{
				if(min(smooth_dx_in_elem,smooth_dy_in_elem)==1)	{
						smooth_elem_patch[e->id]=1;  
						e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("s").marker;
				}
				else{ e->marker = space->get_mesh()->get_element_markers_conversion().get_internal_marker("").marker;
				}
			}
	}



for(int i =0; i<ndof;i++){
			non_smooth = false; smooth_dof[i]=0;
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

	/*Element* elem_neigh=NULL;
	bool p2_neighbor = false;
	for_all_active_elements(e, space->get_mesh()){ 
					for (unsigned int iv = 0; iv < e->get_nvert(); iv++){  
					 	elem_neigh = e->get_neighbor(iv);
						if(elem_neigh!=NULL){ 									
								if(smooth_elem_patch[elem_neigh->id]==0){ 
										if(smooth_elem_patch[e->id]==1){smooth_elem_patch[e->id]=2;break;}
								}
						}
					}
		}*/


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

delete sln; delete R_h_1, delete R_h_2;
}


