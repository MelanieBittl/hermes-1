#ifndef __PREV_SOLUTION_H
#define __PREV_SOLUTION_H
#include "hermes2d.h"
using namespace Hermes;
using namespace Hermes::Hermes2D;

    class PrevSolution : public Solution<double>
    {
		 public:
		   PrevSolution() : Solution<double>(){ };
		    PrevSolution(MeshSharedPtr mesh): Solution<double>(mesh){
							own_mesh = false;     
		    };
		 
		 
		 virtual MeshFunction<double>* clone();      
			void set_own_mesh(MeshSharedPtr mesh);					

		 
	  protected:  
	  bool own_mesh;  


    };
    
class Quad2DCheb;
   
   
 class LimitedSolution : public PrevSolution
 { 
 public:
 LimitedSolution(): PrevSolution() {};

 	LimitedSolution(MeshSharedPtr mesh): PrevSolution(mesh){};

 	void limit_solution(double* coeff_vec, SpaceSharedPtr<double> space); 	

 
 };
    
    
#endif
