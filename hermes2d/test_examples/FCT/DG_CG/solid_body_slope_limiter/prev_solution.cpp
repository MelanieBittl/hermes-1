#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(MeshSharedPtr mesh){						
						if(this->mesh == mesh){
							MeshSharedPtr new_mesh(new Mesh);
							new_mesh->copy(mesh);
							this->mesh = new_mesh;
							own_mesh = true;						
						}else throw Hermes::Exceptions::Exception("Solution mesh unequal own_mesh.");
				}
    
           MeshFunction<double>* PrevSolution::clone(){
                if(this->get_type() == HERMES_SLN)
									return Solution<double>::clone();
								PrevSolution* sln = new PrevSolution(this->mesh);
								return sln;         
          
          };
    
  


//--------------QUAD2dCHEB
  static double3* cheb_tab_tri[11];
    static double3* cheb_tab_quad[11];
    static int      cheb_np_tri[11];
    static int      cheb_np_quad[11];

    static double3** cheb_tab[2] = { cheb_tab_tri, cheb_tab_quad };
    static int*      cheb_np[2]  = { cheb_np_tri,  cheb_np_quad  };

    static class Quad2DCheb : public Quad2D
    {
    public:

      Quad2DCheb()
      {
        max_order[0]  = max_order[1]  = 10;
        num_tables[0] = num_tables[1] = 11;
        tables = cheb_tab;
        np = cheb_np;

        tables[0][0] = tables[1][0] = NULL;
        np[0][0] = np[1][0] = 0;

        int i, j, k, n, m;
        double3* pt;
        for (int mode_i = 0; mode_i <= 1; mode_i++)
        {
          for (k = 0; k <= 10; k++)
          {
            np[mode_i][k] = n = mode_i ? Hermes::sqr(k + 1) : (k + 1)*(k + 2)/2;
            tables[mode_i][k] = pt = new double3[n];

            for (i = k, m = 0; i >= 0; i--)
              for (j = k; j >= (mode_i ? 0 : k-i); j--, m++) {
                pt[m][0] = k ? Hermes::cos(j * M_PI / k) : 1.0;
                pt[m][1] = k ? Hermes::cos(i * M_PI / k) : 1.0;
                pt[m][2] = 1.0;
              }
          }
        }
      };

      ~Quad2DCheb()
      {
        for (int mode_i = 0; mode_i <= 1; mode_i++)
          for (int k = 1; k <= 10; k++)
            delete [] tables[mode_i][k];
      }

      virtual void dummy_fn() {}
    } g_quad_2d_cheb;
    
    
    
        static struct mono_lu_init
    {
    public:

      // this is a set of LU-decomposed matrices shared by all Solutions
      double** mat[2][11];
      int* perm[2][11];

      mono_lu_init()
      {
        memset(mat, 0, sizeof(mat));
      }

      ~mono_lu_init()
      {
        for (int m = 0; m <= 1; m++)
          for (int i = 0; i <= 10; i++)
            if(mat[m][i] != NULL) {
              delete [] mat[m][i];
              delete [] perm[m][i];
            }
      }
    }
    mono_lu;
//--------------------------------------------
  

void LimitedSolution::limit_solution(double* coeff_vec,SpaceSharedPtr<double> space)
{		
		if(space == NULL) throw Exceptions::NullException(1);
      Shapeset *shapeset = space->get_shapeset();
      if(space->get_shapeset() == NULL) throw Exceptions::Exception("Space->shapeset == NULL in limit_solution.");
      PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
      if(pss == NULL) throw Exceptions::Exception("PrecalcShapeset could not be allocated in limit_solution.");
     
  
        if((space->get_shapeset())->get_num_components() == 2) throw Exceptions::Exception("limit_solution not for H_div or H_curl ( num_components() == 2).");
     int o;        
      // Sanity checks.      
      if(space->get_mesh() == NULL) throw Exceptions::Exception("Mesh == NULL in limit_solution.");  
      if(coeff_vec == NULL) throw Exceptions::NullException(3);
      if(coeff_vec == NULL) throw Exceptions::Exception("Coefficient vector == NULL in limit_solution.");
      if(!space->is_up_to_date())
        throw Exceptions::Exception("Provided 'space' is not up to date.");
      if(space->get_shapeset() != pss->get_shapeset())
        throw Exceptions::Exception("Provided 'space' and 'pss' must have the same shapesets.");

			
      free();

      this->space_type = space->get_type();

      this->num_components = pss->get_num_components();
      sln_type = HERMES_SLN;

      // Copy the mesh.
      this->mesh = space->get_mesh();

      // Allocate the coefficient arrays.
      num_elems = this->mesh->get_max_element_id();
      if(elem_orders != NULL)
        delete [] elem_orders;
      elem_orders = new int[num_elems];
      memset(elem_orders, 0, sizeof(int) * num_elems);
      for (int l = 0; l < this->num_components; l++)
      {
        if(elem_coeffs[l] != NULL)
          delete [] elem_coeffs[l];
        elem_coeffs[l] = new int[num_elems];
        memset(elem_coeffs[l], 0, sizeof(int) * num_elems);
      }

      // Obtain element orders, allocate mono_coeffs.
      Element* e;
      num_coeffs = 0;
      for_all_active_elements(e, this->mesh)
      {
        this->mode = e->get_mode();
        o = space->get_element_order(e->id);
        o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
        for (unsigned int k = 0; k < e->get_nvert(); k++)
        {
          int eo = space->get_edge_order(e, k);
          if(eo > o) o = eo;
        }        

        num_coeffs += this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;
        elem_orders[e->id] = o;
      }
      num_coeffs *= this->num_components;
      if(mono_coeffs != NULL)
        delete [] mono_coeffs;
      mono_coeffs = new double[num_coeffs];      
 		int index[4];
      // Express the solution on elements as a linear combination of monomials.
      Quad2D* quad = &g_quad_2d_cheb;
      pss->set_quad_2d(quad);
      double* mono = mono_coeffs;
      for_all_active_elements(e, this->mesh)
      { 
      
		 for (unsigned int iv = 0; iv < e->get_nvert(); iv++)  		
	 			 index[iv] =  shapeset->get_vertex_index(iv,HERMES_MODE_QUAD);
	 			 
        this->mode = e->get_mode();
        o = elem_orders[e->id];
        int np = quad->get_num_points(o, e->get_mode());

        AsmList<double> al;
        space->get_element_assembly_list(e, &al);
        pss->set_active_element(e);

        for (int l = 0; l < this->num_components; l++)
        {
          // Obtain solution values for the current element.
          double* val = mono;
          elem_coeffs[l][e->id] = (int) (mono - mono_coeffs);
          memset(val, 0, sizeof(double)*np);
          for (unsigned int k = 0; k < al.get_cnt(); k++)
          {
            pss->set_active_shape(al.get_idx()[k]);
            pss->set_quad_order(o, H2D_FN_VAL);
            int dof = al.get_dof()[k];
            bool vertex_dof = false;
            for (unsigned int iv = 0; iv < e->get_nvert(); iv++)  
            	if(al.get_idx()[k] == index[iv])
            	{ vertex_dof = true; 
            		break;
            	}
            
            
            double dir_lift_coeff = 0.0;
            // By subtracting space->first_dof we make sure that it does not matter where the
            // enumeration of dofs in the space starts. This ca be either zero or there can be some
            // offset. By adding start_index we move to the desired section of coeff_vec.
            double coef =0.;
            if(vertex_dof == true)
            	coef = al.get_coef()[k] * (dof >= 0 ? coeff_vec[dof] : dir_lift_coeff);
            else
            	coef = al.get_coef()[k] * (dof >= 0 ? 0.0 : dir_lift_coeff);
            			
            double* shape = pss->get_fn_values(l);
            for (int i = 0; i < np; i++)
              val[i] += shape[i] * coef;
          }
          mono += np;

          // solve for the monomial coefficients
          if(mono_lu.mat[this->mode][o] == NULL)
            mono_lu.mat[this->mode][o] = calc_mono_matrix(o, mono_lu.perm[this->mode][o]);
          lubksb(mono_lu.mat[this->mode][o], np, mono_lu.perm[this->mode][o], val);
        }
      }

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL.\n");
      init_dxdy_buffer();
      this->element = NULL;
      if(Solution<double>::static_verbose_output)
        Hermes::Mixins::Loggable::Static::info("Solution: set_coeff_vector - done.");



}


