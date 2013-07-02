#ifndef __H2D_SHAPESET_SERENDIPITY
#define __H2D_SHAPESET_SERENDIPITY


#include "shapeset_h1_all.h" 

namespace Hermes
{
  namespace Hermes2D
  {


   class HERMES_API SerendipityShapeset : public Shapeset
    {
    public:
      SerendipityShapeset();
      virtual Shapeset* clone() { return new SerendipityShapeset (*this); };
      virtual int get_max_index(ElementMode2D mode);
    private:
      virtual int get_id() const { return 4; }
      static const int max_index[H2D_NUM_MODES];
      virtual SpaceType get_space_type() const { return HERMES_L2_SEMI_SPACE; }

      template<typename Scalar> friend class VectorForm;
      template<typename Scalar> friend class MatrixForm;

      template<typename Scalar> friend class DiscreteProblem; 
template<typename Scalar> friend class Solution; 
friend class CurvMap; friend class RefMap; 
template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector; 
template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector; 
template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector; 
template<typename Scalar> friend class RefinementSelectors::OptimumSelector; 
friend class PrecalcShapeset;
    };
}
}


#endif
