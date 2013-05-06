#ifndef __UTIL_H
#define __UTIL_H
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


class VolumeAverage
{
public:
 double calculate_Volume_Average(Element* e,int derivative, Solution<double>* sln_1,Space<double>*ref_space, Solution<double>* sln_2=NULL);
 
 double calculate_norm(Element* e,Space<double>* space, Solution<double>* sln);

	 
};




#endif
