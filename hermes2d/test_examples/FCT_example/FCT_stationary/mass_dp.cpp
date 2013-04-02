#include "mass_dp.h"
    void Mass_DP::assemble_volume_matrix_forms(Stage<double>& stage,
      SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks,
      Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap,
      Hermes::vector<Solution<double>*>& u_ext,
      int marker, Hermes::vector<AsmList<double>*>& al)
    {
      _F_;
      for (unsigned ww = 0; ww < stage.mfvol.size(); ww++)
      {
        MatrixFormVol<double>* mfv = stage.mfvol[ww];
        int m = mfv->i;
        int n = mfv->j;
        if (isempty[m] || isempty[n]) continue;
        if (fabs(mfv->scaling_factor) < 1e-12) continue;

        // Assemble this form only if one of its areas is HERMES_ANY
        // of if the element marker coincides with one of the form's areas.
        bool assemble_this_form = false;
        for (unsigned int ss = 0; ss < mfv->areas.size(); ss++)
        {
          if(mfv->areas[ss] == HERMES_ANY)
          {
            assemble_this_form = true;
            break;
          }
          else
          {
            bool marker_on_space_m = this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_m)
              marker_on_space_m = (this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            bool marker_on_space_n = this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_n)
              marker_on_space_n = (this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            if (marker_on_space_m && marker_on_space_n)
            {
              assemble_this_form = true;
              break;
            }
          }
        }
        if (!assemble_this_form)
          continue;

        // If a block scaling table is provided, and if the scaling get_coef()icient
        // A_mn for this block is zero, then the form does not need to be assembled.
        double block_scaling_coef = 1.;
        if (block_weights != NULL)
        {
          block_scaling_coef = block_weights->get_A(m, n);
          if (fabs(block_scaling_coef) < 1e-12)
            continue;
        }
        bool tra = (m != n) && (mfv->sym != 0);
        bool sym = (m == n) && (mfv->sym == 1);

        // Assemble the local stiffness matrix for the form mfv.
        double **local_stiffness_matrix = NULL;
        local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->get_cnt(), al[n]->get_cnt()));

        for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
        {
          if (!tra && al[m]->get_dof()[i] < 0) continue;
          spss[m]->set_active_shape(al[m]->get_idx()[i]);
          // Unsymmetric block.
          if (!sym)
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all
                  // get_coef()icients multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                  {
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  }
                  local_stiffness_matrix[i][j] = val;
                }
              }
            }
          }
          // Symmetric block.
          else
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              if (j < i && al[n]->get_dof()[j] >= 0)
                continue;
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all get_coef()icients
                  // multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                    refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                }
              }
            }
          }
        }


	//------------- local mass_lumping
			Element* elem_n =refmap[n]->get_active_element();
			Element* elem_m =refmap[m]->get_active_element();
			int elem_id_n =elem_n->id;
			int elem_id_m =elem_m->id;
	if(((spaces[0]->get_element_order(elem_id_m)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_m)==1))&&((spaces[0]->get_element_order(elem_id_n)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_n)==1))){
			Element* elem_neigh =NULL;
			bool p2_neigh = false;
		double elem_diag=elem_n->get_diameter();
			for (unsigned int iv = 0; iv < elem_n->get_nvert(); iv++){  
				elem_neigh = elem_n->get_neighbor(iv);
				if(elem_neigh!=NULL){ 								
					if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
						p2_neigh = true; break;
					}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
				}
				}else if((elem_neigh==NULL)&&(elem_n->en[iv]->bnd==0)){ //kein Randknoten!
						p2_neigh = true; break;
					} 
				if(elem_n->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}

elem_diag=elem_m->get_diameter();
			for (unsigned int iv = 0; iv < elem_m->get_nvert(); iv++){  
				elem_neigh = elem_m->get_neighbor(iv);
					if(elem_neigh!=NULL){ 	info("neigh!=NULL");
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true; break;
						}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
						}
				}else if((elem_neigh==NULL)&&(elem_m->en[iv]->bnd==0)){ //kein Randknoten! 
						p2_neigh = true; break;
				}
				if(elem_m->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}
			if(p2_neigh==false){	//Ordnung soll 1 
					for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
						if(al[m]->get_dof()[i] >= 0)
							for (unsigned int j = 0; j < al[n]->get_cnt(); j++)	
								if(al[n]->get_dof()[j] >= 0)
									if(j!=i){
										 local_stiffness_matrix[i][i] += local_stiffness_matrix[i][j];
										 local_stiffness_matrix[i][j] = 0.;
										}
			}
}



        // Insert the local stiffness matrix into the global one.
        if (mat != NULL)
          mat->add(al[m]->get_cnt(), al[n]->get_cnt(), local_stiffness_matrix, al[m]->get_dof(), al[n]->get_dof());

        // Insert also the off-diagonal (anti-)symmetric block, if required.
        if (tra)
        {
          if (mfv->sym < 0)
            chsgn(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          transpose(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          if (mat != NULL)	
            mat->add(al[n]->get_cnt(), al[m]->get_cnt(), local_stiffness_matrix, al[n]->get_dof(), al[m]->get_dof());


			
        }
      }
    }





//-------------------------Lumped_Flux-----------------------

double max(double a, double b){
			if(a>b) return a;
			else return b;
}

double min(double a, double b){
			if(a<b) return a;
			else return b;
}




    void Lumped_Flux::assemble_volume_matrix_forms(Stage<double>& stage,
      SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks,
      Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap,
      Hermes::vector<Solution<double>*>& u_ext,
      int marker, Hermes::vector<AsmList<double>*>& al)
    {
      _F_;
      for (unsigned ww = 0; ww < stage.mfvol.size(); ww++)
      {
        MatrixFormVol<double>* mfv = stage.mfvol[ww];
        int m = mfv->i;
        int n = mfv->j;
        if (isempty[m] || isempty[n]) continue;
        if (fabs(mfv->scaling_factor) < 1e-12) continue;

        // Assemble this form only if one of its areas is HERMES_ANY
        // of if the element marker coincides with one of the form's areas.
        bool assemble_this_form = false;
        for (unsigned int ss = 0; ss < mfv->areas.size(); ss++)
        {
          if(mfv->areas[ss] == HERMES_ANY)
          {
            assemble_this_form = true;
            break;
          }
          else
          {
            bool marker_on_space_m = this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_m)
              marker_on_space_m = (this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            bool marker_on_space_n = this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_n)
              marker_on_space_n = (this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            if (marker_on_space_m && marker_on_space_n)
            {
              assemble_this_form = true;
              break;
            }
          }
        }
        if (!assemble_this_form)
          continue;

        // If a block scaling table is provided, and if the scaling get_coef()icient
        // A_mn for this block is zero, then the form does not need to be assembled.
        double block_scaling_coef = 1.;
        if (block_weights != NULL)
        {
          block_scaling_coef = block_weights->get_A(m, n);
          if (fabs(block_scaling_coef) < 1e-12)
            continue;
        }
        bool tra = (m != n) && (mfv->sym != 0);
        bool sym = (m == n) && (mfv->sym == 1);

        // Assemble the local stiffness matrix for the form mfv.
        double **local_stiffness_matrix = NULL;
        local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->get_cnt(), al[n]->get_cnt()));

        for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
        {
          if (!tra && al[m]->get_dof()[i] < 0) continue;
          spss[m]->set_active_shape(al[m]->get_idx()[i]);
          // Unsymmetric block.
          if (!sym)
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all
                  // get_coef()icients multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                  {
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  }
                  local_stiffness_matrix[i][j] = val;
                }
              }
            }
          }
          // Symmetric block.
          else
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              if (j < i && al[n]->get_dof()[j] >= 0)
                continue;
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all get_coef()icients
                  // multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                    refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                }
              }
            }
          }
        }


	//------------- local mass_lumping---only for p=1----------
			Element* elem_n =refmap[n]->get_active_element();
			Element* elem_m =refmap[m]->get_active_element();
			int elem_id_n =elem_n->id;
			int elem_id_m =elem_m->id;

	if(((spaces[0]->get_element_order(elem_id_m)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_m)==1))&&((spaces[0]->get_element_order(elem_id_n)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_n)==1))){
			Element* elem_neigh =NULL;
			bool p2_neigh = false;
		double elem_diag=elem_n->get_diameter();

			for (unsigned int iv = 0; iv < elem_n->get_nvert(); iv++){  
				elem_neigh = elem_n->get_neighbor(iv);
				if(elem_neigh!=NULL){ 								
					if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
						p2_neigh = true; break;
					}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
				}
				}else if((elem_neigh==NULL)&&(elem_n->en[iv]->bnd==0)){ //kein Randknoten!
						p2_neigh = true; break;
					} 
				if(elem_n->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}

elem_diag=elem_m->get_diameter();
			for (unsigned int iv = 0; iv < elem_m->get_nvert(); iv++){  
				elem_neigh = elem_m->get_neighbor(iv);
					if(elem_neigh!=NULL){ 	info("neigh!=NULL");
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true; break;
						}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
						}
				}else if((elem_neigh==NULL)&&(elem_m->en[iv]->bnd==0)){ //kein Randknoten! 
						p2_neigh = true; break;
				}
				if(elem_m->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}

			if(p2_neigh==false){	//Ordnung soll 1 				
				for (int i = 0; i < al[m]->get_cnt(); i++){						
						if(i>=4) printf("mass_dp: i:444444444444444");
						if(al[m]->get_dof()[i] >= 0){
							elem_flux[i] = 0;
							elem_flux[i] -=local_stiffness_matrix[i][i]*u_H[al[m]->get_dof()[i]];
							for (int j = 0; j < al[n]->get_cnt(); j++)	{
								if(j>=4) printf("mass_dp: i:444444444444444");
								if(al[n]->get_dof()[j] >= 0)
									if(j!=i){
										 local_stiffness_matrix[i][i] += local_stiffness_matrix[i][j];
										elem_flux[i]-=local_stiffness_matrix[i][j]*u_H[al[n]->get_dof()[j]];  //-M_c u_h
										 local_stiffness_matrix[i][j] = 0.;
										}
							}
							elem_flux[i]+=local_stiffness_matrix[i][i]*u_H[al[m]->get_dof()[i]]; //f = (M_L - M_c)u_h
						P_plus[al[m]->get_dof()[i]] += max(elem_flux[i],0.);
						P_minus[al[m]->get_dof()[i]] += min(elem_flux[i],0.);
						for (int j = 0; j < al[n]->get_cnt(); j++){
Q_plus[al[m]->get_dof()[i]] = max(Q_plus[al[m]->get_dof()[i]], local_stiffness_matrix[i][i]*(u_L[al[n]->get_dof()[j]]-u_L[al[m]->get_dof()[i]]));
Q_minus[al[m]->get_dof()[i]] = min(Q_minus[al[m]->get_dof()[i]], local_stiffness_matrix[i][i]*(u_L[al[n]->get_dof()[j]]-u_L[al[m]->get_dof()[i]]));
}
						}
				}
		}
}
//---------------------------------------------------
        // Insert the local stiffness matrix into the global one.
   /*     if (mat != NULL)
          mat->add(al[m]->get_cnt(), al[n]->get_cnt(), local_stiffness_matrix, al[m]->get_dof(), al[n]->get_dof());

        // Insert also the off-diagonal (anti-)symmetric block, if required.
        if (tra)
        {
          if (mfv->sym < 0)
            chsgn(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          transpose(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          if (mat != NULL)	
            mat->add(al[n]->get_cnt(), al[m]->get_cnt(), local_stiffness_matrix, al[n]->get_dof(), al[m]->get_dof());


			
        }*/
      }
    }


   void Lumped_Flux::assemble_volume_vector_forms(Stage<double>& stage,
      SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks,
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss,
      Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution<double>*>& u_ext,
      int marker, Hermes::vector<AsmList<double>*>& al)
    {
      _F_;

      if(rhs == NULL)
        return;
			double alpha_elem =1.;
for (unsigned int i = 0; i < 4; i++) elem_flux[i] =0.;

  for (unsigned ww = 0; ww < stage.mfvol.size(); ww++)
      {
        MatrixFormVol<double>* mfv = stage.mfvol[ww];
        int m = mfv->i;
        int n = mfv->j;
        if (isempty[m] || isempty[n]) continue;
        if (fabs(mfv->scaling_factor) < 1e-12) continue;

        // Assemble this form only if one of its areas is HERMES_ANY
        // of if the element marker coincides with one of the form's areas.
        bool assemble_this_form = false;
        for (unsigned int ss = 0; ss < mfv->areas.size(); ss++)
        {
          if(mfv->areas[ss] == HERMES_ANY)
          {
            assemble_this_form = true;
            break;
          }
          else
          {
            bool marker_on_space_m = this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_m)
              marker_on_space_m = (this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            bool marker_on_space_n = this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_n)
              marker_on_space_n = (this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            if (marker_on_space_m && marker_on_space_n)
            {
              assemble_this_form = true;
              break;
            }
          }
        }
        if (!assemble_this_form)
          continue;

        // If a block scaling table is provided, and if the scaling get_coef()icient
        // A_mn for this block is zero, then the form does not need to be assembled.
        double block_scaling_coef = 1.;
        if (block_weights != NULL)
        {
          block_scaling_coef = block_weights->get_A(m, n);
          if (fabs(block_scaling_coef) < 1e-12)
            continue;
        }
        bool tra = (m != n) && (mfv->sym != 0);
        bool sym = (m == n) && (mfv->sym == 1);

        // Assemble the local stiffness matrix for the form mfv.
        double **local_stiffness_matrix = NULL;
        local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->get_cnt(), al[n]->get_cnt()));

        for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
        {
          if (!tra && al[m]->get_dof()[i] < 0) continue;
          spss[m]->set_active_shape(al[m]->get_idx()[i]);
          // Unsymmetric block.
          if (!sym)
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all
                  // get_coef()icients multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                  {
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  }
                  local_stiffness_matrix[i][j] = val;

                }
              }
            }
          }
          // Symmetric block.
          else
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              if (j < i && al[n]->get_dof()[j] >= 0)
                continue;
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all get_coef()icients
                  // multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                    refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                }
              }
            }
          }
        }


	//------------- local mass_lumping---only for p=1----------
			Element* elem_n =refmap[n]->get_active_element();
			Element* elem_m =refmap[m]->get_active_element();
			int elem_id_n =elem_n->id;
			int elem_id_m =elem_m->id;

	if(((spaces[0]->get_element_order(elem_id_m)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_m)==1))&&((spaces[0]->get_element_order(elem_id_n)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_n)==1))){
			Element* elem_neigh =NULL;
			bool p2_neigh = false;

		double elem_diag=elem_n->get_diameter();

			for (unsigned int iv = 0; iv < elem_n->get_nvert(); iv++){  
				elem_neigh = elem_n->get_neighbor(iv);
				if(elem_neigh!=NULL){ 								
					if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
						p2_neigh = true; break;
					}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
				}
				}else if((elem_neigh==NULL)&&(elem_n->en[iv]->bnd==0)){ //kein Randknoten!
						/*elem_neigh= elem_n->parent->get_neighbor(iv);
						if(elem_neigh!=NULL)
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true; break;
						}*/
						p2_neigh = true; break;
					} 
				if(elem_n->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}

elem_diag=elem_m->get_diameter();
			for (unsigned int iv = 0; iv < elem_m->get_nvert(); iv++){  
				elem_neigh = elem_m->get_neighbor(iv);
					if(elem_neigh!=NULL){ 	info("neigh!=NULL");
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true; break;
						}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
						}
				}else if((elem_neigh==NULL)&&(elem_m->en[iv]->bnd==0)){ //kein Randknoten! 
						/*elem_neigh= elem_m->parent->get_neighbor(iv);
						if(elem_neigh!=NULL)
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true;  break;
						}*/
						p2_neigh = true; break;
				}
				if(elem_m->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}
			if(p2_neigh==false){	//Ordnung soll 1 		
				for (int i = 0; i < al[m]->get_cnt(); i++){	
						double loc =	local_stiffness_matrix[i][i];
						//elem_flux[i]=0.;			
						if(i>=4) printf("mass_dp: i:444444444444444");		
						if(al[m]->get_dof()[i] >= 0){
							elem_flux[i] =-local_stiffness_matrix[i][i]*u_H[al[m]->get_dof()[i]];
							for (int j = 0; j < al[n]->get_cnt(); j++)	{
								if(j>=4) printf("mass_dp: i:444444444444444");
								if(al[n]->get_dof()[j] >= 0)
									if(al[m]->get_dof()[i]!=al[n]->get_dof()[i])  printf("mass_dp: dofs not equal");
									if(j!=i){
										 local_stiffness_matrix[i][i] = local_stiffness_matrix[i][i]+local_stiffness_matrix[i][j];									
											elem_flux[i]=elem_flux[i]-local_stiffness_matrix[i][j]*u_H[al[n]->get_dof()[j]];  //-M_c u_h
										 local_stiffness_matrix[i][j] = 0.;
										}
									
							}			
					

							elem_flux[i]=elem_flux[i]+local_stiffness_matrix[i][i]*u_H[al[m]->get_dof()[i]]; //f = (M_L - M_c)u_h

						}
				}
				for (int i = 0; i < al[m]->get_cnt(); i++){
						if(al[m]->get_dof()[i] >= 0){
										if(elem_flux[i]>0) alpha_elem = R_plus[al[m]->get_dof()[i]];
										else alpha_elem = R_minus[al[m]->get_dof()[i]];
										for (int j = 0; j < al[n]->get_cnt(); j++){
												if(al[n]->get_dof()[j] >= 0){
																if(elem_flux[j]>0) alpha_elem = min(alpha_elem,R_plus[al[n]->get_dof()[j]]);
																else alpha_elem = min(alpha_elem,R_minus[al[n]->get_dof()[j]]);
												}
										}
					}
				} 


				
		}
}
        // Insert the local stiffness matrix into the global one.
       if (mat != NULL)
          mat->add(al[m]->get_cnt(), al[n]->get_cnt(), local_stiffness_matrix, al[m]->get_dof(), al[n]->get_dof());

        // Insert also the off-diagonal (anti-)symmetric block, if required.
        if (tra)
        {
          if (mfv->sym < 0)
            chsgn(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          transpose(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          if (mat != NULL)	
            mat->add(al[n]->get_cnt(), al[m]->get_cnt(), local_stiffness_matrix, al[n]->get_dof(), al[m]->get_dof());


			
        }
      }

//---------------------------------------------------
      for (unsigned int ww = 0; ww < stage.vfvol.size(); ww++)
      {
        VectorFormVol<double>* vfv = stage.vfvol[ww];
        int m = vfv->i;
        if (isempty[vfv->i]) continue;
        if (fabs(vfv->scaling_factor) < 1e-12) continue;

        // Assemble this form only if one of its areas is HERMES_ANY
        // of if the element marker coincides with one of the form's areas.
        bool assemble_this_form = false;
        for (unsigned int ss = 0; ss < vfv->areas.size(); ss++)
        {
          if(vfv->areas[ss] == HERMES_ANY)
          {
            assemble_this_form = true;
            break;
          }
          else
          {
            bool marker_on_space_m = this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(vfv->areas[ss]).valid;
            if(marker_on_space_m)
              marker_on_space_m = (this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(vfv->areas[ss]).marker == marker);

            if (marker_on_space_m)
            {
              assemble_this_form = true;
              break;
            }
          }
        }
        if (assemble_this_form == false) continue;

        double* vector = NULL;
        
        for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
        {
          if (al[m]->get_dof()[i] < 0) continue;

          spss[m]->set_active_shape(al[m]->get_idx()[i]);

          // Numerical integration performed only if the coefficient
          // multiplying the form is nonzero.
					//add flux correction
          if (std::abs(al[m]->get_coef()[i]) > 1e-12){
            rhs->add(al[m]->get_dof()[i], (eval_form(vfv, u_ext, spss[m], refmap[m])+alpha_elem*elem_flux[i]) * al[m]->get_coef()[i]);
					}
        }
      }

    }




 void Lumped_Flux::assemble(double* coeff_vec, SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks, bool add_dir_lift,
      Table* block_weights)
    {
      _F_

      // Check that the block scaling table have proper dimension.
      if (block_weights != NULL)
        if (block_weights->get_size() != wf->get_neq())
          throw Exceptions::LengthException(6, block_weights->get_size(), wf->get_neq());

      // Creating matrix sparse structure.
      create_sparse_structure(mat, rhs, force_diagonal_blocks, block_weights);

      // Convert the coefficient vector into vector of external solutions.
      Hermes::vector<Solution<double>*> u_ext;
      int first_dof = 0;
      if (coeff_vec != NULL) for (int i = 0; i < wf->get_neq(); i++)
      {
        Solution<double>* external_solution_i = new Solution<double>(spaces[i]->get_mesh());
        Solution<double>::vector_to_solution(coeff_vec, spaces[i], external_solution_i, add_dir_lift, first_dof);
        u_ext.push_back(external_solution_i);
        first_dof += spaces[i]->get_num_dofs();
      }

      // Reset the warnings about insufficiently high integration order.
      reset_warn_order();

      // Create slave pss's, refmaps.
      Hermes::vector<PrecalcShapeset *> spss;
      Hermes::vector<RefMap *> refmap;

      // Initialize slave pss's, refmaps.
      initialize_psss(spss);
      initialize_refmaps(refmap);

      // Initialize matrix buffer.
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
      if (mat != NULL)
        get_matrix_buffer(9);

      // Create assembling stages.
      Hermes::vector<Stage<double> > stages = Hermes::vector<Stage<double> >();
      bool want_matrix = (mat != NULL);
      bool want_vector = (rhs != NULL);
      wf->get_stages(spaces, u_ext, stages, want_matrix, want_vector);

      // Loop through all assembling stages -- the purpose of this is increased performance
      // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
      // In such a case, the matrix forms are assembled over one mesh, and only the rhs
      // traverses through the union mesh. On the other hand, if you don't use multi-mesh
      // at all, there will always be only one stage in which all forms are assembled as usual.
      for (unsigned ss = 0; ss < stages.size(); ss++)
        // Assemble one stage. One stage is a collection of functions,
        // and meshes that can not be further minimized.
        // E.g. if a linear form uses two external solutions, each of
        // which is defined on a different mesh, and different to the
        // mesh of the current test function, then the stage would have
        // three meshes. By stage functions, all functions are meant: shape
        // functions (their precalculated values), and mesh functions.
        assemble_one_stage(stages[ss], mat, NULL, force_diagonal_blocks,
        block_weights, spss, refmap, u_ext);

			for(int i =0;i<ref_ndof; i++){
					if(P_plus[i]!=0.) R_plus[i] = min(1.,Q_plus[i]/P_plus[i]) ;
					else 							 R_plus[i] = 1.;
					if(P_minus[i]!=0.) R_minus[i] = min(1.,Q_minus[i]/P_minus[i]);
					else 							 R_minus[i] = 1.;
			}


      // Loop through all assembling stages -- the purpose of this is increased performance
      // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
      // In such a case, the matrix forms are assembled over one mesh, and only the rhs
      // traverses through the union mesh. On the other hand, if you don't use multi-mesh
      // at all, there will always be only one stage in which all forms are assembled as usual.
      for (unsigned ss = 0; ss < stages.size(); ss++)
        // Assemble one stage. One stage is a collection of functions,
        // and meshes that can not be further minimized.
        // E.g. if a linear form uses two external solutions, each of
        // which is defined on a different mesh, and different to the
        // mesh of the current test function, then the stage would have
        // three meshes. By stage functions, all functions are meant: shape
        // functions (their precalculated values), and mesh functions.
        assemble_one_stage(stages[ss], mat, rhs, force_diagonal_blocks,
        block_weights, spss, refmap, u_ext);



      // Deinitialize matrix buffer.
      if(matrix_buffer != NULL)
        delete [] matrix_buffer;
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;

      // Deinitialize slave pss's, refmaps.
      for(Hermes::vector<PrecalcShapeset *>::iterator it = spss.begin(); it != spss.end(); it++)
        delete *it;
      for(Hermes::vector<RefMap *>::iterator it = refmap.begin(); it != refmap.end(); it++)
        delete *it;

      // Delete the vector u_ext.
      for(typename Hermes::vector<Solution<double>*>::iterator it = u_ext.begin(); it != u_ext.end(); it++)
        delete *it;
    }


    void Lumped_Flux::assemble(double* coeff_vec, Vector<double>* rhs,
      bool force_diagonal_blocks, bool add_dir_lift,
      Table* block_weights)
    {
      _F_;
      assemble(coeff_vec, NULL, rhs, force_diagonal_blocks, add_dir_lift, block_weights);
    }

    void Lumped_Flux::assemble(SparseMatrix<double>* mat, Vector<double>* rhs,
      bool force_diagonal_blocks, Table* block_weights)
    {
      _F_;
      double* coeff_vec = NULL;
      assemble(coeff_vec, mat, rhs, force_diagonal_blocks, block_weights);
    }

//----------Artificial Diffusion-------------------



 void Artificial_Diffusion ::assemble_volume_matrix_forms(Stage<double>& stage,
      SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks,
      Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap,
      Hermes::vector<Solution<double>*>& u_ext,
      int marker, Hermes::vector<AsmList<double>*>& al)
    {
      _F_;
      for (unsigned ww = 0; ww < stage.mfvol.size(); ww++)
      {
        MatrixFormVol<double>* mfv = stage.mfvol[ww];
        int m = mfv->i;
        int n = mfv->j;
        if (isempty[m] || isempty[n]) continue;
        if (fabs(mfv->scaling_factor) < 1e-12) continue;

        // Assemble this form only if one of its areas is HERMES_ANY
        // of if the element marker coincides with one of the form's areas.
        bool assemble_this_form = false;
        for (unsigned int ss = 0; ss < mfv->areas.size(); ss++)
        {
          if(mfv->areas[ss] == HERMES_ANY)
          {
            assemble_this_form = true;
            break;
          }
          else
          {
            bool marker_on_space_m = this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_m)
              marker_on_space_m = (this->spaces[m]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            bool marker_on_space_n = this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).valid;
            if(marker_on_space_n)
              marker_on_space_n = (this->spaces[n]->get_mesh()->get_element_markers_conversion().get_internal_marker(mfv->areas[ss]).marker == marker);

            if (marker_on_space_m && marker_on_space_n)
            {
              assemble_this_form = true;
              break;
            }
          }
        }
        if (!assemble_this_form)
          continue;

        // If a block scaling table is provided, and if the scaling get_coef()icient
        // A_mn for this block is zero, then the form does not need to be assembled.
        double block_scaling_coef = 1.;
        if (block_weights != NULL)
        {
          block_scaling_coef = block_weights->get_A(m, n);
          if (fabs(block_scaling_coef) < 1e-12)
            continue;
        }
        bool tra = (m != n) && (mfv->sym != 0);
        bool sym = (m == n) && (mfv->sym == 1);

        // Assemble the local stiffness matrix for the form mfv.
        double **local_stiffness_matrix = NULL;
        local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->get_cnt(), al[n]->get_cnt()));

        for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
        {
          if (!tra && al[m]->get_dof()[i] < 0) continue;
          spss[m]->set_active_shape(al[m]->get_idx()[i]);
          // Unsymmetric block.
          if (!sym)
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all
                  // get_coef()icients multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                  {
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  }
                  local_stiffness_matrix[i][j] = val;
                }
              }
            }
          }
          // Symmetric block.
          else
          {
            for (unsigned int j = 0; j < al[n]->get_cnt(); j++)
            {
              if (j < i && al[n]->get_dof()[j] >= 0)
                continue;
              pss[n]->set_active_shape(al[n]->get_idx()[j]);
              if (al[n]->get_dof()[j] >= 0)
              {
                if (mat != NULL)
                {
                  double val = 0;
                  // Numerical integration performed only if all get_coef()icients
                  // multiplying the form are nonzero.
                  if (std::abs(al[m]->get_coef()[i]) > 1e-12 && std::abs(al[n]->get_coef()[j]) > 1e-12)
                    val = block_scaling_coef * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                    refmap[m]) * al[n]->get_coef()[j] * al[m]->get_coef()[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                }
              }
            }
          }
        }


	//------------- local mass_lumping
			Element* elem_n =refmap[n]->get_active_element();
			Element* elem_m =refmap[m]->get_active_element();
			int elem_id_n =elem_n->id;
			int elem_id_m =elem_m->id;
	if(((spaces[0]->get_element_order(elem_id_m)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_m)==1))&&((spaces[0]->get_element_order(elem_id_n)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_n)==1))){
			Element* elem_neigh =NULL;
			bool p2_neigh = false;
		double elem_diag=elem_n->get_diameter();
			for (unsigned int iv = 0; iv < elem_n->get_nvert(); iv++){  
				elem_neigh = elem_n->get_neighbor(iv);
				if(elem_neigh!=NULL){ 								
					if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
						p2_neigh = true; break;
					}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
				}
				}else if((elem_neigh==NULL)&&(elem_n->en[iv]->bnd==0)){ //kein Randknoten!
						p2_neigh = true; break;
					} 
				if(elem_n->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}

elem_diag=elem_m->get_diameter();
			for (unsigned int iv = 0; iv < elem_m->get_nvert(); iv++){  
				elem_neigh = elem_m->get_neighbor(iv);
					if(elem_neigh!=NULL){ 	info("neigh!=NULL");
						if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
							p2_neigh = true; break;
						}else if(elem_diag != elem_neigh->get_diameter()){	// Nachbar ist kleiner => haengender Knoten, kein FCT ueber haengendne Knoten hinweg
											p2_neigh =true; 
											break;
						}
				}else if((elem_neigh==NULL)&&(elem_m->en[iv]->bnd==0)){ //kein Randknoten! 
						p2_neigh = true; break;
				}
				if(elem_m->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
			}
			if(p2_neigh==false){	//Ordnung soll 1 
			double a=0.; double b=0.;;
					for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
						if(al[m]->get_dof()[i] >= 0){
									local_stiffness_matrix[i][i]=0;
							for (unsigned int j = (i+1); j < al[n]->get_cnt(); j++)	
								if(al[n]->get_dof()[j] >= 0)
									if(j!=i){
														a= -local_stiffness_matrix[i][j];
														b=-local_stiffness_matrix[j][i];
												if((a>=b)&&(a>0.0)){
															local_stiffness_matrix[i][j] = a;
															local_stiffness_matrix[j][i] = a;
															local_stiffness_matrix[i][i] -=a;
												}else if((b>a)&&(b>0.0)){
															local_stiffness_matrix[i][j] = b;
															local_stiffness_matrix[j][i] = b;
															local_stiffness_matrix[i][i] -=b;
												}else{
															local_stiffness_matrix[i][j] = 0.;
															local_stiffness_matrix[j][i] = 0.;
												}

										}
					}
			}
}



        // Insert the local stiffness matrix into the global one.
        if (mat != NULL)
          mat->add(al[m]->get_cnt(), al[n]->get_cnt(), local_stiffness_matrix, al[m]->get_dof(), al[n]->get_dof());

        // Insert also the off-diagonal (anti-)symmetric block, if required.
        if (tra)
        {
          if (mfv->sym < 0)
            chsgn(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          transpose(local_stiffness_matrix, al[m]->get_cnt(), al[n]->get_cnt());

          if (mat != NULL)	
            mat->add(al[n]->get_cnt(), al[m]->get_cnt(), local_stiffness_matrix, al[n]->get_dof(), al[m]->get_dof());


			
        }
      }
    }


