#include "prev_solution.h"
      

		void PrevSolution::set_own_mesh(const Mesh* mesh)
		{						
				if(this->mesh == mesh){
					Mesh* new_mesh = new Mesh;
					new_mesh->copy(mesh);
					this->mesh = new_mesh;
					own_mesh = true;						
				}else throw Hermes::Exceptions::Exception("Solution mesh unequal own_mesh.");
		}

	  MeshFunction<double>* PrevSolution::clone()
	  {
						PrevSolution* sln = new PrevSolution;						
						 sln->copy(this);
						// sln->set_own_mesh(this->mesh);
						return sln;         
	 
	 };
    

  
    void PrevSolution::free()
    {
        
		if((own_mesh == true)&&(mesh!=NULL)){ own_mesh = false; delete mesh; }
		mesh = NULL;
		Solution<double>::free();

    }




  
