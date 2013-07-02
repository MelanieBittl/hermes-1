
#include "shapeset_serendipity.h" 
#include "global.h"
#include "shapeset.h"
#include "shapeset_common.h"

#define s_1(x) (1-x)
#define s_2(x) (1+x)
#define s_3(x,y) (x-y-1)
#define s_4(x,y) (-x-y-1)
#define s_5(x,y) (x+y-1)
#define s_6(x) (1-(x*x))



//// quad serendipity shapeset /////////////////////////////////////////////////////////////////
namespace Hermes
{
  namespace Hermes2D
  {

    static double serendipity_1(double x, double y)
    {
      return  0.25*s_1(x)*s_1(y)*s_4(x,y);
    }
		    static double serendipity_1x(double x, double y)
    {
      return   0.25*s_1(y)*(2*x+y);
    }
	    static double serendipity_1y(double x, double y)
    {
      return  0.25*s_1(x)*(2*y+x);
    }
    static double serendipity_1xx(double x, double y)
    {
      return  0.5*s_1(y) ;
    }
    static double serendipity_1xy(double x, double y)
    {
      return 0.25*(1-2*x-2*y);
    }
    static double serendipity_1yy(double x, double y)
    {
      return 0.5*s_1(x);
    }


    static double serendipity_2(double x, double y)
    {
      return  0.25*s_2(x)*s_1(y)*s_3(x,y);
    }
		    static double serendipity_2x(double x, double y)
    {
      return   0.25*s_1(y)*(2*x-y);
    }
	    static double serendipity_2y(double x, double y)
    {
      return  0.25*s_2(x)*(2*y-x);
    }
    static double serendipity_2xx(double x, double y)
    {
      return  0.5*s_1(y) ;
    }
    static double serendipity_2xy(double x, double y)
    {
      return 0.25*(-1-2*x+2*y);
    }
    static double serendipity_2yy(double x, double y)
    {
      return 0.5*s_2(x);
    }


    static double serendipity_3(double x, double y)
    {
      return  0.25*s_2(x)*s_2(y)*s_5(x,y);
    }
		    static double serendipity_3x(double x, double y)
    {
      return   0.25*s_2(y)*(2*x+y);
    }
	    static double serendipity_3y(double x, double y)
    {
      return  0.25*s_2(x)*(2*y+x);
    }
    static double serendipity_3xx(double x, double y)
    {
      return  0.5*s_2(y) ;
    }
    static double serendipity_3xy(double x, double y)
    {
      return 0.25*(1+2*x+2*y);
    }
    static double serendipity_3yy(double x, double y)
    {
      return 0.5*s_2(x);
    }


    static double serendipity_4(double x, double y)
    {
      return  0.25*s_1(x)*s_2(y)*s_3(y,x);
    }
		    static double serendipity_4x(double x, double y)
    {
      return   0.25*s_2(y)*(2*x-y);
    }
	    static double serendipity_4y(double x, double y)
    {
      return  0.25*s_1(x)*(2*y-x);
    }
    static double serendipity_4xx(double x, double y)
    {
      return  0.5*s_2(y) ;
    }
    static double serendipity_4xy(double x, double y)
    {
      return 0.25*(-1+2*x-2*y);
    }
    static double serendipity_4yy(double x, double y)
    {
      return 0.5*s_1(x);
    }



    static double serendipity_5(double x, double y)
    {
      return  0.5*s_6(x)*s_1(y);
    }
		    static double serendipity_5x(double x, double y)
    {
      return   -x*s_1(y);
    }
	    static double serendipity_5y(double x, double y)
    {
      return  -0.5*s_6(x);
    }
    static double serendipity_5xx(double x, double y)
    {
      return -s_1(y);
    }
    static double serendipity_5xy(double x, double y)
    {
      return x;
    }
    static double serendipity_5yy(double x, double y)
    {
      return 0;
    }




    static double serendipity_6(double x, double y)
    {
      return  0.5*s_6(y)*s_2(x);
    }
		    static double serendipity_6x(double x, double y)
    {
      return   0.5*s_6(y);
    }
	    static double serendipity_6y(double x, double y)
    {
      return  0.5*s_2(x)*(-2*y);
    }
    static double serendipity_6xx(double x, double y)
    {
      return 0;
    }
    static double serendipity_6xy(double x, double y)
    {
      return -y;
    }
    static double serendipity_6yy(double x, double y)
    {
      return 0.5*s_2(x)*(-2*y);
    }


    static double serendipity_7(double x, double y)
    {
      return  0.5*s_6(x)*s_2(y);
    }
		    static double serendipity_7x(double x, double y)
    {
      return   s_2(y)*(-x);
    }
	    static double serendipity_7y(double x, double y)
    {
      return  0.5*s_6(x);
    }
    static double serendipity_7xx(double x, double y)
    {
      return -s_2(y);
    }
    static double serendipity_7xy(double x, double y)
    {
      return -x;
    }
    static double serendipity_7yy(double x, double y)
    {
      return 0;
    }


    static double serendipity_8(double x, double y)
    {
      return  0.5*s_6(y)*s_1(x);
    }
		    static double serendipity_8x(double x, double y)
    {
      return   0.5*(y*y-1);
    }
	    static double serendipity_8y(double x, double y)
    {
      return  -y*s_1(x);
    }
    static double serendipity_8xx(double x, double y)
    {
      return 0;
    }
    static double serendipity_8xy(double x, double y)
    {
      return y;
    }
    static double serendipity_8yy(double x, double y)
    {
      return -s_1(x);
    }


    static Shapeset::shape_fn_t serendipity_quad_fn[] =
    {
      serendipity_1,   serendipity_2,   serendipity_3, serendipity_4, serendipity_5, serendipity_6, serendipity_7, serendipity_8 
    };
    static Shapeset::shape_fn_t serendipity_quad_fn_dx[] =
    {
      serendipity_1x,   serendipity_2x,   serendipity_3x, serendipity_4x, serendipity_5x, serendipity_6x, serendipity_7x, serendipity_8x 
    };
    static Shapeset::shape_fn_t serendipity_quad_fn_dy[] =
    {
      serendipity_1y,   serendipity_2y,   serendipity_3y, serendipity_4y, serendipity_5y, serendipity_6y, serendipity_7y, serendipity_8y 
    };
    static Shapeset::shape_fn_t serendipity_quad_fn_dxx[] =
    {
      serendipity_1xx,   serendipity_2xx,   serendipity_3xx, serendipity_4xx, serendipity_5xx, serendipity_6xx, serendipity_7xx, serendipity_8xx 
    };
    static Shapeset::shape_fn_t serendipity_quad_fn_dxy[] =
    {
      serendipity_1xy,   serendipity_2xy,   serendipity_3xy, serendipity_4xy, serendipity_5xy, serendipity_6xy, serendipity_7xy, serendipity_8xy 
    };
    static Shapeset::shape_fn_t serendipity_quad_fn_dyy[] =
    {
      serendipity_1yy,   serendipity_2yy,   serendipity_3yy, serendipity_4yy, serendipity_5yy, serendipity_6yy, serendipity_7yy, serendipity_8yy 
    };

    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table[1]     = { serendipity_quad_fn };
    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table_dx[1]  = { serendipity_quad_fn_dx };
    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table_dy[1]  = { serendipity_quad_fn_dy };
    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table_dxx[1] = { serendipity_quad_fn_dxx };
    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table_dxy[1] = { serendipity_quad_fn_dxy };
    Shapeset::shape_fn_t* serendipity_quad_shape_fn_table_dyy[1] = { serendipity_quad_fn_dyy };

  
    #define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL
   int* serendipity_quad_bubble_indices[] =
    {
      NULL16, NULL16,
      NULL16, NULL16,
      NULL16, NULL16,
      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,
      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16
    };
  #define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    int serendipity_quad_bubble_count[] =
    {
      zero16, zero16,
      zero16, zero16,
     zero16, zero16,
		zero16,zero16,		zero16,zero16,		zero16,zero16,		zero16,zero16,
		zero16,zero16,		zero16,zero16,		zero16,zero16,		zero16,zero16
    };

    int serendipity_quad_vertex_indices[4] = { 0, 1, 2, 3 };

    static int serendipity_quad_edge_indices_0[6] =  {  0, 1, 1, 0, 4, 4 };
    static int serendipity_quad_edge_indices_1[6] =  { 	1, 2, 2, 1, 5, 5 };
    static int serendipity_quad_edge_indices_2[6] =  { 	2, 3, 3, 2, 6, 6 };
    static int serendipity_quad_edge_indices_3[6] =  {  3, 0, 0, 3, 7, 7 };

    int* serendipity_quad_edge_indices[4] =
    {
      serendipity_quad_edge_indices_0,
      serendipity_quad_edge_indices_1,
      serendipity_quad_edge_indices_2,
      serendipity_quad_edge_indices_3
    };

#define oo H2D_MAKE_QUAD_ORDER
#define XX(a, b) oo(a, b), oo(a, b)

    int serendipity_quad_index_to_order[] =
    {
      oo(2, 2),   oo(2, 2),   oo(2, 2),  
      oo(2, 2),   oo(2, 1),   oo(1, 2),  
      oo(2, 1),   oo(1, 2)
     
    };



    static Shapeset::shape_fn_t** serendipity_shape_fn_table[2] =
    {
      NULL,
      serendipity_quad_shape_fn_table
    };

    static Shapeset::shape_fn_t** serendipity_shape_fn_table_dx[2] =
    {
      NULL,
      serendipity_quad_shape_fn_table_dx
    };

    static Shapeset::shape_fn_t** serendipity_shape_fn_table_dy[2] =
    {
      NULL,
      serendipity_quad_shape_fn_table_dy
    };

    static Shapeset::shape_fn_t** serendipity_shape_fn_table_dxx[2] =
    {
      NULL,
      serendipity_quad_shape_fn_table_dxx
    };

    static Shapeset::shape_fn_t** serendipity_shape_fn_table_dyy[2] =
    {
     NULL,
      serendipity_quad_shape_fn_table_dyy
    };

    static Shapeset::shape_fn_t** serendipity_shape_fn_table_dxy[2] =
    {
      NULL,
      serendipity_quad_shape_fn_table_dxy
    };

    static int* serendipity_vertex_indices[2] =
    {
      NULL,
      serendipity_quad_vertex_indices
    };

    static int** serendipity_edge_indices[2] =
    {
      NULL,
      serendipity_quad_edge_indices
    };

    static int** serendipity_bubble_indices[2] =
    {
      NULL,
      serendipity_quad_bubble_indices
    };

    static int* serendipity_bubble_count[2] =
    {
      NULL,
      serendipity_quad_bubble_count
    };

    static int* serendipity_index_to_order[2] =
    {
      NULL,
      serendipity_quad_index_to_order
    };
      SerendipityShapeset::SerendipityShapeset()
		{
      shape_table[0] = serendipity_shape_fn_table;
      shape_table[1] = serendipity_shape_fn_table_dx;
      shape_table[2] = serendipity_shape_fn_table_dy;
      shape_table[3] = serendipity_shape_fn_table_dxx;
      shape_table[4] = serendipity_shape_fn_table_dyy;
      shape_table[5] = serendipity_shape_fn_table_dxy;

      vertex_indices = serendipity_vertex_indices;
      edge_indices = serendipity_edge_indices;
      bubble_indices = serendipity_bubble_indices;
      bubble_count = serendipity_bubble_count;
      index_to_order = serendipity_index_to_order;

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
      min_order = 2;
      num_components = 1;

      ebias = 2;

      comb_table = NULL;

		}

    int SerendipityShapeset::get_max_index(ElementMode2D mode) { return max_index[mode]; }

    const int  SerendipityShapeset::max_index[2] = { 0, 7 };
}


}

