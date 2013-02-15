#include "prev_solution.h"
      
      
				void PrevSolution::set_own_mesh(const Mesh* mesh){						
						if(this->mesh == mesh){
							Mesh* new_mesh = new Mesh;
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
    
 void PrevSolution::copy(const Solution<double>* sln)
 {
	 if((own_mesh == true)&&(mesh!=NULL)) { delete mesh; own_mesh = false; mesh =NULL;}
 				
 	 Solution<double>::copy(sln);
 
 }
  
    void PrevSolution::free()
    {        
		if((own_mesh == true)&&(mesh!=NULL)){ delete mesh; own_mesh = false; mesh =NULL;}
		mesh = NULL; 
 				
      Solution<double>::free();
    }


