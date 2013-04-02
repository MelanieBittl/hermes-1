#ifndef __HP_ADAPT_H
#define __HP_ADAPT_H
#include "hermes2d.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


class  HPAdapt : public Adapt<double>
{
public:
     HPAdapt(Hermes::vector<Space<double>*> spaces, int max_p,double h_min,double h_max, int max_dof):Adapt<double>(spaces)
     {
    this->h_min = h_min;
     this->max_p= max_p; this->h_max=h_max; this->max_dof=max_dof;
     
};
	  HPAdapt(Space<double>* space, ProjNormType proj_norm): Adapt<double>(space, proj_norm){	
	} ; 
	~HPAdapt(){};
 
bool adapt(int* elements_to_refine,int* no_of_refinement_steps,  int* smooth_elem=NULL, int regularize = -1);


protected:

double h_min;
double h_max;
int max_p;
int max_dof;

};







#endif
