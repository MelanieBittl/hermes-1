#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
     HPAdapt(Hermes::vector<SpaceSharedPtr<double> > spaces,ErrorCalculator<double>* error_calculator):Adapt<double>(spaces,error_calculator){
};
	  HPAdapt(SpaceSharedPtr<double> space, ErrorCalculator<double>* error_calculator, AdaptivityStoppingCriterion<double>* strategy =NULL): Adapt<double>(space, error_calculator, strategy){	
	} ; 
	~HPAdapt(){};
 
bool adapt(int* elements_to_refine,int* no_of_refinement_steps, double h_min,double h_max, int max_dof, int regularize = -1);


bool reduce_order(std::set<int>* discont_elem);



};







#endif
