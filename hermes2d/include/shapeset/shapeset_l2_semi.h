#ifndef __H2D_SHAPESET_L2_SEMI_ALL
#define __H2D_SHAPESET_L2_SEMI_ALL


#include "shapeset_h1_all.h"
namespace Hermes
{
  namespace Hermes2D
  {


    class HERMES_API L2SEMIShapeset : public H1ShapesetJacobi
    {
    public:
      virtual Shapeset* clone() { return new L2SEMIShapeset (*this); };
    private:
      virtual int get_id() const { return 3; }
			virtual SpaceType get_space_type() const { return HERMES_L2_SEMI_SPACE; }

      template<typename Scalar> friend class VectorForm;
      template<typename Scalar> friend class MatrixForm;      

      template<typename Scalar> friend class DiscreteProblem; 
			template<typename Scalar> friend class Solution; 
			friend class CurvMap; 
			friend class RefMap; 
			template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector; 
			template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector; 
			template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector; 
			template<typename Scalar> friend class RefinementSelectors::OptimumSelector; 
			friend class PrecalcShapeset;
    };

  }
}
#endif
