#ifndef __H2D_SHAPESET_BB
#define __H2D_SHAPESET_BB


#include "shapeset.h"
namespace Hermes
{
  namespace Hermes2D
  {


    /// @ingroup spaces
    /// Shape functions based on Bernstein-Bezier polynoms.
    class HERMES_API ShapesetBB : public Shapeset
    {
    public:
      ShapesetBB(int order);
    typedef double (*shape_fn_bb)(double, double, int);
   virtual double get_value(int n, int index, double x, double y, int component, ElementMode2D mode);

      virtual Shapeset* clone() { return new ShapesetBB(space_order); };
      virtual SpaceType get_space_type() const { return HERMES_H1_SPACE; }
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 5; }

		protected:
			int space_order;

			shape_fn_bb*** shape_table_bb[6];
      /// Constructs the linear combination of edge functions, forming a constrained edge function.
      ///
      virtual double get_constrained_value(int n, int index, double x, double y, int component, ElementMode2D mode);
      
      static const int max_index[H2D_NUM_MODES];
    };

  }
}
#endif
