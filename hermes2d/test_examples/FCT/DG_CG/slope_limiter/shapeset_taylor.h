#ifndef __SHAPESET_TAYLOR_H
#define __SHAPESET_TAYLOR_H
#include "shapeset.h"
#include "hermes2d.h" 
#include "util.h" 

using namespace Hermes;
using namespace Hermes::Hermes2D;

class AuxiliarySolution : public ExactSolutionScalar<double>
{
public:
  AuxiliarySolution(MeshSharedPtr mesh, bool x_component, double c) : ExactSolutionScalar<double>(mesh), x_component(x_component), c(c) {};
 

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;

   MeshFunction<double>* clone() const ;
   bool x_component;
   double c;
   
};


//Shapeset of local Taylor basis functions

 class HERMES_API TaylorShapeset : public Shapeset
 {
 public:
    typedef double (*shape_fn_t)(double, double, double, double , double);
    TaylorShapeset();
    TaylorShapeset(Element* e,SpaceSharedPtr<double>  space);
    ~TaylorShapeset();
   void init_data(Element* e,SpaceSharedPtr<double>  new_space);
        
   virtual Shapeset* clone() { return new TaylorShapeset(active_element,space); };
   
   virtual SpaceType get_space_type() const { return HERMES_L2_SPACE; }
   virtual int get_max_index(ElementMode2D mode);
   
   virtual double get_value(int n, int index, double x, double y, int component, ElementMode2D mode);
    
    double get_delta_x(){return delta_x;};
    double get_delta_y(){return delta_y;};
  
   
 protected:
   virtual int get_id() const { return 32; }

   static const int max_index[2];
   
   Element* active_element;
   SpaceSharedPtr<double> space;
   MeshSharedPtr mesh;
   double x_c;
   double y_c;
   double x_max, x_min, y_max, y_min;
   double delta_x, delta_y;
   VolumeAverage* vol_av;
   AuxiliarySolution* sln_x;
   AuxiliarySolution* sln_y;
   double rest_taylor_4;
   double rest_taylor_5;
   double rest_taylor_6;

      shape_fn_t*** shape_table[6];

   template<typename Scalar> friend class DiscreteProblem;
   template<typename Scalar> friend class VectorForm;
   template<typename Scalar> friend class MatrixForm;
   template<typename Scalar> friend class Solution;
   friend class CurvMap; friend class RefMap;
   template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector;
   template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector;
   template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector;
   template<typename Scalar> friend class RefinementSelectors::OptimumSelector; friend class PrecalcShapeset;




 };







#endif
