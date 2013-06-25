#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
	  HPAdapt(SpaceSharedPtr<double> space, ErrorCalculator<double>* error_calculator, AdaptivityStoppingCriterion<double>* strategy =NULL): Adapt<double>(space, error_calculator, strategy){	
	} ; 
	~HPAdapt(){};
 
bool adapt(int* elements_to_refine,int* no_of_refinement_steps, int max_p,double h_min,double h_max, int max_dof, int regularize = -1);

bool adapt_smooth(int* smooth_dof, int max_p);
};







#endif
