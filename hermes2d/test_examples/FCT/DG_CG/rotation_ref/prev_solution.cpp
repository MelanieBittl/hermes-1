#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(MeshSharedPtr mesh){						
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
								PrevSolution* sln = new PrevSolution();
								sln->copy(this);
								return sln;    
          
          };
    

  
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




  
