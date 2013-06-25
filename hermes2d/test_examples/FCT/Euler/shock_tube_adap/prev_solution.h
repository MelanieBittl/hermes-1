#ifndef __PREV_SOLUTION_H
#define __PREV_SOLUTION_H
#include "hermes2d.h"
using namespace Hermes;
using namespace Hermes::Hermes2D;

    class PrevSolution : public Solution<double>
    {
    public:
      PrevSolution() : Solution<double>()
      {
      	this->mesh = NULL;
      	own_mesh =false; 
      };

    
    ~PrevSolution()
    {
      if((own_mesh == true)&&(mesh!=NULL)) { own_mesh = false; delete mesh; }		
      mesh = NULL;
    }
    
    virtual MeshFunction<double>* clone();      
		void set_own_mesh(const Mesh* mesh);					

    
    
  protected:  
  bool own_mesh;  
    void free();




    };
    
    
    
#endif
