
#include "global.h"
#include "shapeset_common.h"
#include "shapeset_bb.h"

namespace Hermes
{
  namespace Hermes2D
  {
    // ORDER 1

    // Vertex functions, order 1

    // number 1
    inline double bb_f1(double x, double y, int order=0)
    {
      return std::pow(lambda2(x, y),order);
    }

    inline double bb_f1_dx(double x, double y,int order=0)
    {
      return std::pow(lambda2(x, y), order-1)* lambda2x(x,y)*order;
    }

    inline double bb_f1_dy(double x, double y, int order=0)
    {
      return std::pow(lambda2(x, y), order-1)*lambda2y(x, y)*order;
    }

    inline double bb_f1_dxx(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda2(x, y), order-2)*lambda2x(x, y)*lambda2x(x, y)));
			else return 0.;
    }

    inline double bb_f1_dyy(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda2(x, y), order-2)*lambda2y(x, y)*lambda2y(x, y)));
			else return 0.;
    }

    inline double bb_f1_dxy(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda2(x, y), order-2)*lambda2x(x, y)*lambda2y(x, y)));
			else return 0.;
    }

    // number 2
    inline double bb_f2(double x, double y, int order=0)
    {
      return std::pow(lambda3(x, y),order);
    }

    inline double bb_f2_dx(double x, double y, int order=0)
    {
      return std::pow(lambda3(x, y), order-1)* lambda3x(x,y)*order;
    }

    inline double bb_f2_dy(double x, double y, int order=0)
    {
      return std::pow(lambda3(x, y), order-1)* lambda3y(x,y)*order;
    }

    inline double bb_f2_dxx(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda3(x, y), order-2)*lambda3x(x, y)*lambda3x(x, y)));
			else return 0.;
    }

    inline double bb_f2_dyy(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda3(x, y), order-2)*lambda3y(x, y)*lambda3y(x, y)));
			else return 0.;
    }

    inline double bb_f2_dxy(double x, double y, int order=0)
    {
			if(order>1)
      	return (order*((order-1)*std::pow(lambda3(x, y), order-2)*lambda3x(x, y)*lambda3y(x, y)));
			else return 0.;
    }

    // number 3
    inline double bb_f3(double x, double y, int order=0)
    {
        return std::pow(lambda1(x, y),order);
    }

    inline double bb_f3_dx(double x, double y, int order=0)
    {
      return std::pow(lambda1(x, y), order-1)* lambda1x(x,y)*order;
    }

    inline double bb_f3_dy(double x, double y, int order=0)
    {
      return std::pow(lambda1(x, y), order-1)* lambda1y(x,y)*order;
    }

    inline double bb_f3_dxx(double x, double y, int order=0)
    {
     	if(order>1)
      	return (order*((order-1)*std::pow(lambda1(x, y), order-2)*lambda1x(x, y)*lambda1x(x, y)));
			else return 0.;
    }

    inline double bb_f3_dyy(double x, double y, int order=0)
    {
     	if(order>1)
      	return (order*((order-1)*std::pow(lambda1(x, y), order-2)*lambda1y(x, y)*lambda1y(x, y)));
			else return 0.;
    }

    inline double bb_f3_dxy(double x, double y, int order=0)
    {
     			if(order>1)
      	return (order*((order-1)*std::pow(lambda1(x, y), order-2)*lambda1y(x, y)*lambda1x(x, y)));
			else return 0.;
    }

    // ORDER 2

    // Edge functions, order 2

    // number 4
    inline double bb_f4(double x, double y, int order=0)
    {
				return (2.*bb_f1(x,y,order-1)*bb_f2(x,y,order-1));
    }

    inline double bb_f4_dx(double x, double y, int order=0)
    {
			return (2.*(bb_f1_dx(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dx(x,y,order-1)*bb_f1(x,y,order-1)));
    }

    inline double bb_f4_dy(double x, double y, int order=0)
    {
				return (2.*(bb_f1_dy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dy(x,y,order-1)*bb_f1(x,y,order-1)));
    }

    inline double bb_f4_dxx(double x, double y, int order=0)
    {
			return (2.*(bb_f1_dxx(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dxx(x,y,order-1)*bb_f1(x,y,order-1)+2*bb_f1_dx(x,y,order-1)*bb_f2_dx(x,y,order-1)));
    }

    inline double bb_f4_dyy(double x, double y, int order=0)
    {
			return (2.*(bb_f1_dyy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dyy(x,y,order-1)*bb_f1(x,y,order-1)+2*bb_f1_dy(x,y,order-1)*bb_f2_dy(x,y,order-1)));
    }

    inline double bb_f4_dxy(double x, double y, int order=0)
    {
			return (2.*(bb_f1_dxy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dxy(x,y,order-1)*bb_f1(x,y,order-1)+bb_f1_dx(x,y,order-1)*bb_f2_dy(x,y,order-1)+bb_f1_dy(x,y,order-1)*bb_f2_dx(x,y,order-1)));
    }

    // number 5
    inline double bb_f5(double x, double y, int order=0)
    {
				return (2.*bb_f2(x,y,order-1)*bb_f3(x,y,order-1));
    }

    inline double bb_f5_dx(double x, double y, int order=0)
    {
			return 2.*(bb_f2_dx(x,y,order-1)*bb_f3(x,y,order-1)+bb_f3_dx(x,y,order-1)*bb_f2(x,y,order-1));
    }

    inline double bb_f5_dy(double x, double y, int order=0)
    {
				return 2.*(bb_f3_dy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dy(x,y,order-1)*bb_f3(x,y,order-1));
    }

    inline double bb_f5_dxx(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dxx(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dxx(x,y,order-1)*bb_f3(x,y,order-1)+2*bb_f3_dx(x,y,order-1)*bb_f2_dx(x,y,order-1));
    }

    inline double bb_f5_dyy(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dyy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dyy(x,y,order-1)*bb_f3(x,y,order-1)+2*bb_f3_dy(x,y,order-1)*bb_f2_dy(x,y,order-1));
    }

    inline double bb_f5_dxy(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dxy(x,y,order-1)*bb_f2(x,y,order-1)+bb_f2_dxy(x,y,order-1)*bb_f3(x,y,order-1)+bb_f3_dx(x,y,order-1)*bb_f2_dy(x,y,order-1)+bb_f3_dy(x,y,order-1)*bb_f2_dx(x,y,order-1));
    }

    // number 6
    inline double bb_f6(double x, double y, int order=0)
    {
				return (2.*bb_f1(x,y,order-1)*bb_f3(x,y,order-1));
    }

    inline double bb_f6_dx(double x, double y, int order=0)
    {
			return 2.*(bb_f1_dx(x,y,order-1)*bb_f3(x,y,order-1)+bb_f3_dx(x,y,order-1)*bb_f1(x,y,order-1));
    }

    inline double bb_f6_dy(double x, double y, int order=0)
    {
				return 2.*(bb_f3_dy(x,y,order-1)*bb_f1(x,y,order-1)+bb_f1_dy(x,y,order-1)*bb_f3(x,y,order-1));
    }

    inline double bb_f6_dxx(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dxx(x,y,order-1)*bb_f1(x,y,order-1)+bb_f1_dxx(x,y,order-1)*bb_f3(x,y,order-1)+2*bb_f3_dx(x,y,order-1)*bb_f1_dx(x,y,order-1));
    }

    inline double bb_f6_dyy(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dyy(x,y,order-1)*bb_f1(x,y,order-1)+bb_f1_dyy(x,y,order-1)*bb_f3(x,y,order-1)+2*bb_f3_dy(x,y,order-1)*bb_f1_dy(x,y,order-1));
    }

    inline double bb_f6_dxy(double x, double y, int order=0)
    {
			return 2.*(bb_f3_dxy(x,y,order-1)*bb_f1(x,y,order-1)+bb_f1_dxy(x,y,order-1)*bb_f3(x,y,order-1)+bb_f3_dx(x,y,order-1)*bb_f1_dy(x,y,order-1)+bb_f3_dy(x,y,order-1)*bb_f1_dx(x,y,order-1));
    }

 

    ////////////////////////////////////////////////////////////////////////////////////////////////////

    static ShapesetBB::shape_fn_bb bb_tri_fn[] =
    {
      bb_f1,    bb_f2,    bb_f3,    bb_f4,    bb_f5,    bb_f6
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dx[] =
    {
      bb_f1_dx,   bb_f2_dx,   bb_f3_dx,   bb_f4_dx,   bb_f5_dx,   bb_f6_dx
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dy[] =
    {
      bb_f1_dy,   bb_f2_dy,   bb_f3_dy,   bb_f4_dy,   bb_f5_dy,   bb_f6_dy
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dxx[] =
    {
      bb_f1_dxx,   bb_f2_dxx,   bb_f3_dxx,   bb_f4_dxx,   bb_f5_dxx,   bb_f6_dxx
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dyy[] =
    {
      bb_f1_dyy,   bb_f2_dyy,   bb_f3_dyy,   bb_f4_dyy,   bb_f5_dyy,   bb_f6_dyy
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dxy[] =
    {
      bb_f1_dxy,   bb_f2_dxy,   bb_f3_dxy,   bb_f4_dxy,   bb_f5_dxy,   bb_f6_dxy
    };


    static int* bb_tri_bubble_indices[3] =
    {
      NULL, NULL, NULL,
    };

    static int bb_tri_bubble_count[3] = { 0, 0, 0};

    static int bb_tri_edge_indices_0[6] =  { 0, 1, 1, 0, 3, 3 };
    static int bb_tri_edge_indices_1[6] =  { 1, 2, 2, 1, 4, 4};
    static int bb_tri_edge_indices_2[6] =  { 2, 0, 0, 2, 5, 5 };

    static int* bb_tri_edge_indices[3] =
    {
      bb_tri_edge_indices_0,
      bb_tri_edge_indices_1,
      bb_tri_edge_indices_2
    };

    static int bb_tri_vertex_indices[3] = { 0, 1, 2 };

    static int bb_tri_index_to_order[6] =
    {
      1, 1, 1, 2, 2, 2
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table[1] =
    {
      bb_tri_fn
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table_dx[1] =
    {
      bb_tri_fn_dx
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table_dy[1] =
    {
      bb_tri_fn_dy
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table_dxx[1] =
    {
      bb_tri_fn_dxx
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table_dyy[1] =
    {
      bb_tri_fn_dyy
    };

    static ShapesetBB::shape_fn_bb* bb_tri_shape_fn_table_dxy[1] =
    {
      bb_tri_fn_dxy
    };

////////////////////////

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table[2] =
    {
      bb_tri_shape_fn_table,
NULL
    };

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table_dx[2] =
    {
      bb_tri_shape_fn_table_dx,
      NULL
    };

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table_dy[2] =
    {
      bb_tri_shape_fn_table_dy,
      NULL
    };

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table_dxx[2] =
    {
      bb_tri_shape_fn_table_dxx,
      NULL
    };

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table_dyy[2] =
    {
      bb_tri_shape_fn_table_dyy,
      NULL
    };

    static ShapesetBB::shape_fn_bb** bb_shape_fn_table_dxy[2] =
    {
      bb_tri_shape_fn_table_dxy,
      NULL
    };

    static int* bb_vertex_indices[2] =
    {
      bb_tri_vertex_indices,
     NULL
    };

    static int** bb_edge_indices[2] =
    {
      bb_tri_edge_indices,
     NULL
    };

    static int** bb_bubble_indices[2] =
    {
      bb_tri_bubble_indices,
      NULL
    };

    static int* bb_bubble_count[2] =
    {
      bb_tri_bubble_count,
      NULL
    };

    static int* bb_index_to_order[2] =
    {
      bb_tri_index_to_order,
      NULL
    };

    int ShapesetBB::get_max_index(ElementMode2D mode) { return max_index[mode]; }

    ShapesetBB::ShapesetBB(int order)
    {
			this->space_order = order;

      shape_table_bb[0] = bb_shape_fn_table;
      shape_table_bb[1] = bb_shape_fn_table_dx;
      shape_table_bb[2] = bb_shape_fn_table_dy;
      shape_table_bb[3] = bb_shape_fn_table_dxx;
      shape_table_bb[4] = bb_shape_fn_table_dyy;
      shape_table_bb[5] = bb_shape_fn_table_dxy;

      vertex_indices = bb_vertex_indices;
      edge_indices = bb_edge_indices;
      bubble_indices = bb_bubble_indices;
      bubble_count = bb_bubble_count;
      index_to_order = bb_index_to_order;

      ref_vert[0][0][0] = -1.0;
      ref_vert[0][0][1] = -1.0;
      ref_vert[0][1][0] =  1.0;
      ref_vert[0][1][1] = -1.0;
      ref_vert[0][2][0] = -1.0;
      ref_vert[0][2][1] =  1.0;

      ref_vert[1][0][0] = -1.0;
      ref_vert[1][0][1] = -1.0;
      ref_vert[1][1][0] =  1.0;
      ref_vert[1][1][1] = -1.0;
      ref_vert[1][2][0] =  1.0;
      ref_vert[1][2][1] =  1.0;
      ref_vert[1][3][0] = -1.0;
      ref_vert[1][3][1] =  1.0;

      max_order = 2;
      min_order = 1;
      num_components = 1;

      ebias = 2;

      comb_table = NULL;
    }

double ShapesetBB::get_value(int n, int index, double x, double y, int component, ElementMode2D mode)
{     
	 if(index >= 0)
      {
        ShapesetBB::shape_fn_bb** shape_expansion = shape_table_bb[n][mode];
        if(shape_expansion == NULL)
        { // requested expansion (f, df/dx, df/dy, ddf/dxdx, ...) is not defined.
          //just to keep the number of warnings low: warn just once about a given combinations of n, mode, and index.
          static int warned_mode = -1, warned_index = -1, warned_n = 1;
          this->warn_if(warned_mode != mode || warned_index != index || warned_n != n, "Requested undefined expansion %d (mode: %d) of a shape %d, returning 0", n, mode, index);
          warned_mode = mode;
          warned_index = index;
          warned_n = n;
          return 0.;
        }
        else
        {  
      
         return shape_expansion[component][index](x, y, space_order);
       
        }
      }
      else
        return get_constrained_value(n, index, x, y, component, mode);
}
    double ShapesetBB::get_constrained_value(int n, int index, double x, double y, int component, ElementMode2D mode)
    {
      index = -1 - index;

      int part = (unsigned) index >> 7;
      int order = (index >> 3) & 15;
      int edge = (index >> 1) & 3;
      int ori = index & 1;

      int i, nc;
      double sum, *comb = get_constrained_edge_combination(order, part, ori, nc, mode);

      sum = 0.0;
      shape_fn_bb* table = shape_table_bb[n][mode][component];
      for (i = 0; i < nc; i++)
        sum += comb[i] * table[get_edge_index(edge, ori, i + ebias, mode)](x, y,space_order);

      return sum;
    }



    const int ShapesetBB::max_index[2] = { 5, 0 };
  }
}
