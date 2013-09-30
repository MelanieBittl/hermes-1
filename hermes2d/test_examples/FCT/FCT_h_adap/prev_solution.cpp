#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(const MeshSharedPtr mesh){						
						if(this->mesh == mesh)
            {
							MeshSharedPtr new_mesh(new Mesh);
							new_mesh->copy(mesh);
							this->mesh = new_mesh;
							own_mesh = true;						
						}
            else
              throw Hermes::Exceptions::Exception("Solution mesh unequal own_mesh.");
				}
    
           MeshFunction<double>* PrevSolution::clone(){
                if(this->get_type() == HERMES_SLN)
									return Solution<double>::clone();
								PrevSolution* sln = new PrevSolution(this->mesh);
								return sln;         
          
          };
    
  void PrevSolution::copy(const MeshFunction<double>* sln)
	{

      const Solution<double>* solution = dynamic_cast<const Solution<double>*>(sln);
      if(solution == NULL)
        throw Exceptions::Exception("The instance is in fact not a Solution instance in copy().");

      if(solution->get_type() == HERMES_UNDEF)
        throw Hermes::Exceptions::Exception("Solution being copied is uninitialized.");
      free();

			if(own_mesh) this->mesh->copy(solution->get_mesh());
      else {
			this->mesh = solution->get_mesh();
			this->set_own_mesh(solution->get_mesh());
			}

      sln_type = solution->get_type();
      space_type = solution->get_space_type();
      this->num_components = solution->get_num_components();
      num_dofs = solution->get_num_dofs();

     if(solution->get_type() == HERMES_SLN) // standard solution: copy coefficient arrays
      {
        num_coeffs = solution->get_num_coeffs();
        num_elems = solution->get_num_elems();

        mono_coeffs = new double[num_coeffs];
        memcpy(mono_coeffs, solution->mono_coeffs, sizeof(double) * num_coeffs);

        for (int l = 0; l < this->num_components; l++)
        {
          elem_coeffs[l] = new int[num_elems];
          memcpy(elem_coeffs[l], solution->get_elem_coeffs(l), sizeof(int) * num_elems);
        }

        elem_orders = new int[num_elems];
        memcpy(elem_orders, solution->elem_orders, sizeof(int) * num_elems);

        init_dxdy_buffer();
      }
      else // Const, exact handled differently.
        throw Hermes::Exceptions::Exception("Undefined or exact solutions cannot be copied into an instance of Solution already coming from computation.");

      this->element = NULL;
	}


  
    void PrevSolution::free()
    {
      if(mono_coeffs  != NULL) { delete [] mono_coeffs;   mono_coeffs = NULL;  }
      if(elem_orders != NULL) { delete [] elem_orders;  elem_orders = NULL; }
      if(dxdy_buffer != NULL) { delete [] dxdy_buffer;  dxdy_buffer = NULL; }

      for (int i = 0; i < this->num_components; i++)
        if(elem_coeffs[i] != NULL)
        { delete [] elem_coeffs[i];  elem_coeffs[i] = NULL; }
        
 				
        e_last = NULL;

        free_tables();
    }




  
