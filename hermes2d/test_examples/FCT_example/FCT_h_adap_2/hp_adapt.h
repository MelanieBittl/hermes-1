#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
	  HPAdapt(Space<double>* space, ProjNormType proj_norm): Adapt<double>(space, proj_norm){	
	} ; 
	~HPAdapt(){};
 
bool adapt(int* elements_to_refine,int* no_of_refinement_steps, int max_p,double h_min,double h_max, int max_dof, int regularize = -1);


};







#endif
