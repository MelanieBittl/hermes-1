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


	// local mass_lumping
	Element* elem_n =refmap[n]->get_active_element();
	Element* elem_m =refmap[m]->get_active_element();
	int elem_id_n =elem_n->id;
	int elem_id_m =elem_m->id;
	Element* elem_neigh =NULL;
	bool p2_neigh = false;

	for (unsigned int iv = 0; iv < elem_n->get_nvert(); iv++){  
		elem_neigh = elem_n->get_neighbor(iv);
		if(elem_neigh!=NULL) 
			if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
				p2_neigh = true; break;
			}
		if(elem_n->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
	}
		for (unsigned int iv = 0; iv < elem_m->get_nvert(); iv++){  
		elem_neigh = elem_m->get_neighbor(iv);
			if(elem_neigh!=NULL) 
				if((spaces[0]->get_element_order(elem_neigh->id)== H2D_MAKE_QUAD_ORDER(2, 2))||(spaces[0]->get_element_order(elem_neigh->id)==2)){
					p2_neigh = true; break;
				}
		if(elem_m->vn[iv]->is_constrained_vertex() ==true){p2_neigh = true; break;}
	}
	if(((spaces[0]->get_element_order(elem_id_m)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_m)==1))&&((spaces[0]->get_element_order(elem_id_n)== H2D_MAKE_QUAD_ORDER(1, 1))||(spaces[0]->get_element_order(elem_id_n)==1))&&(p2_neigh==false)){	//Ordnung soll 1 
	    for (unsigned int i = 0; i < al[m]->get_cnt(); i++)
				if(al[m]->get_dof()[i] >= 0)
				  for (unsigned int j = 0; j < al[n]->get_cnt(); j++)	
						if(al[n]->get_dof()[j] >= 0)
							if(j!=i){
								 local_stiffness_matrix[i][i] += local_stiffness_matrix[i][j];
								 local_stiffness_matrix[i][j] = 0.;
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

