#ifndef __PREV_SOLUTION_H
#define __PREV_SOLUTION_H
#include "hermes2d.h"
using namespace Hermes;
using namespace Hermes::Hermes2D;

    class PrevSolution : public Solution<double>
    {
		 public:
		   PrevSolution() : Solution<double>(){ };
		    PrevSolution(const Mesh* mesh): Solution<double>(mesh){
							own_mesh = false;     
		    };
		 
		 ~PrevSolution()
		 {
		   if((own_mesh == true)&&(mesh!=NULL)) { delete mesh; own_mesh = false; mesh =NULL;}
		   mesh = NULL;
		 }
		 
		 virtual MeshFunction<double>* clone();      
			void set_own_mesh(const Mesh* mesh);					

		  virtual void copy(const Solution<double>* sln);
		 
	  protected:  
	  bool own_mesh;  
	  virtual void free();

    };

    
#endif
