#ifndef __H2D_SHAPESET_BB
#define __H2D_SHAPESET_BB


// Common definitions used by the shapesets...
#define bc1(x,y) (-((x) + (y)) / 2.)
#define bc2(x,y) (((x) + 1.) / 2.)
#define bc3(x,y) (((y) + 1.) / 2.)

/// x derivatives of affine coordinates
#define bc1x(x,y) (-1.0 / 2.0)
#define bc2x(x,y) (1.0 / 2.0)
#define bc3x(x,y) (0.0)

/// y derivatives of affine coordinates
#define bc1y(x,y) (-1.0 / 2.0)
#define bc2y(x,y) (0.0)
#define bc3y(x,y) (1.0 / 2.0)



// Common definitions used by the shapesets...
#define lam1(x) ((1.-x)/2.)
#define lam2(x) (((x) + 1.) / 2.)


/// x derivatives of affine coordinates
#define lam1d(x) (-1.0 / 2.0)
#define lam2d(x) (1.0 / 2.0)






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
			~ShapesetBB();
    typedef double (*shape_fn_bb)(double, double, int);
   virtual double get_value(int n, int index, double x, double y, int component, ElementMode2D mode);

virtual double get_fn_value_order(int order, int index, double x, double y, int component, ElementMode2D mode);

      virtual Shapeset* clone() { return new ShapesetBB(space_order); };
      virtual SpaceType get_space_type() const { return HERMES_H1_SPACE; }
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 5; }

		protected:
			int space_order;
			//int no_shape_fn;
			//int no_shape_fn_quad;

			shape_fn_bb*** shape_table_bb[6];
      /// Constructs the linear combination of edge functions, forming a constrained edge function.
      ///
      virtual double get_constrained_value(int n, int index, double x, double y, int component, ElementMode2D mode);
      
    };

  }
}
#endif
