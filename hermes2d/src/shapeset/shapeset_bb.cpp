
#include "global.h"
#include "shapeset_common.h"
#include "shapeset_bb.h"

namespace Hermes
{
  namespace Hermes2D
  {

		int faculty(int i)
		{
				int result=1.;
				while(i>1){ result *=i; i--;}
				return result;
		}

		inline double B(int order,int i, int j, int k, double x, double y)
		{			
				if(order<0) return 0.;
				else{
					double frac = faculty(order)/(faculty(i)*faculty(j)*faculty(k));
					return frac*(std::pow(bc1(x,y),i)*std::pow(bc2(x,y),j)*std::pow(bc3(x,y),k));	
				}	
		}

		inline double Bx(int order, int i, int j, int k, double x, double y)
		{		if(order>0)		
				{	double ix = 0; double jx = 0; double kx = 0;
					if(i>0) ix =	order*B(order-1,i-1, j, k, x, y)*bc1x(x,y);
					if(j>0) jx = 	order*B(order-1,i, j-1, k, x, y)*bc2x(x,y);
					if(k>0) kx = 	order*B(order-1,i, j, k-1, x, y)*bc3x(x,y);
					return (ix+jx+kx);	
				}	else return 0.;
		}

		inline double By(int order, int i, int j, int k, double x, double y)
		{
				if(order>0)	{			
					double iy = 0; double jy = 0; double ky = 0;
					if(i>0) iy =	order*B(order-1,i-1, j, k, x, y)*bc1y(x,y);
					if(j>0) jy = 	order*B(order-1,i, j-1, k, x, y)*bc2y(x,y);
					if(k>0) ky = 	order*B(order-1,i, j, k-1, x, y)*bc3y(x,y);
					return (iy+jy+ky);	
				}else return 0.;	
		}


		inline double Bxx(int order, int i, int j, int k, double x, double y)
		{		
				if(order>1)	{		
					double ix = 0; double jx = 0; double kx = 0;
					if(i>0) ix =	bc1x(x,y)*Bx(order-1,i-1,j,k,x,y) ;
					if(j>0) jx = 	bc2x(x,y)*Bx(order-1,i,j-1,k,x,y);
					if(k>0) kx = 	bc3x(x,y)*Bx(order-1,i,j,k-1,x,y);
					return order*(ix+jx+kx);	
				}else return 0.;		
		}

		inline double Byy(int order, int i, int j, int k, double x, double y)
		{
				if(order>1)	{		
					double iy = 0; double jy = 0; double ky = 0;
					if(i>0) iy =	bc1y(x,y)*By(order-1,i-1,j,k,x,y) ;
					if(j>0) jy = 	bc2y(x,y)*By(order-1,i,j-1,k,x,y);
					if(k>0) ky = 	bc3y(x,y)*By(order-1,i,j,k-1,x,y);
					return order*(iy+jy+ky);	
				}else return 0.;	
		}	

		inline double Bxy(int order, int i, int j, int k, double x, double y)
		{				
				if(order>1)	{		
					double ix = 0; double jx = 0; double kx = 0;
					if(i>0) ix =	bc1x(x,y)*By(order-1,i-1,j,k,x,y) ;
					if(j>0) jx = 	bc2x(x,y)*By(order-1,i,j-1,k,x,y);
					if(k>0) kx = 	bc3x(x,y)*By(order-1,i,j,k-1,x,y);
					return order*(ix+jx+kx);
				}else return 0.;			
		}






    // ORDER 1
    // Vertex functions

    // number 1
    inline double bb_f1(double x, double y, int order=0)
    {
			return B(order,order,0,0,x,y);
    }

    inline double bb_f1_dx(double x, double y,int order=0)
    {
      return Bx(order,order,0,0,x,y);
    }

    inline double bb_f1_dy(double x, double y, int order=0)
    {
      return By(order,order,0,0,x,y);
    }

    inline double bb_f1_dxx(double x, double y, int order=0)
    {
      return Bxx(order,order,0,0,x,y);
    }

    inline double bb_f1_dyy(double x, double y, int order=0)
    {
      return Byy(order,order,0,0,x,y);
    }

    inline double bb_f1_dxy(double x, double y, int order=0)
    {
      return Bxy(order,order,0,0,x,y);
    }

    // number 2
    inline double bb_f2(double x, double y, int order=0)
    {
      return B(order,0,order,0,x,y);
    }

    inline double bb_f2_dx(double x, double y, int order=0)
    {
     return Bx(order,0,order,0,x,y);
    }

    inline double bb_f2_dy(double x, double y, int order=0)
    {
      return By(order,0,order,0,x,y);
    }

    inline double bb_f2_dxx(double x, double y, int order=0)
    {
			return Bxx(order,0,order,0,x,y);
    }

    inline double bb_f2_dyy(double x, double y, int order=0)
    {
      return Byy(order,0,order,0,x,y);
    }

    inline double bb_f2_dxy(double x, double y, int order=0)
    {
      return Bxy(order,0,order,0,x,y);
    }

    // number 3
    inline double bb_f3(double x, double y, int order=0)
    {
      return B(order,0,0,order,x,y);
    }

    inline double bb_f3_dx(double x, double y, int order=0)
    {
      return Bx(order,0,0,order,x,y);
    }

    inline double bb_f3_dy(double x, double y, int order=0)
    {
      return By(order,0,0,order,x,y);
    }

    inline double bb_f3_dxx(double x, double y, int order=0)
    {
      return Bxx(order,0,0,order,x,y);
    }

    inline double bb_f3_dyy(double x, double y, int order=0)
    {
      return Byy(order,0,0,order,x,y);
    }

    inline double bb_f3_dxy(double x, double y, int order=0)
    {
      return Bxy(order,0,0,order,x,y);
    }

    // ORDER 2
    // Edge functions

    // number 4
    inline double bb_f4(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return B(order,i,j,0,x,y);
    }

    inline double bb_f4_dx(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bx(order,i,j,0,x,y);
    }

    inline double bb_f4_dy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return By(order,i,j,0,x,y);
    }

    inline double bb_f4_dxx(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bxx(order,i,j,0,x,y);
    }

    inline double bb_f4_dyy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Byy(order,i,j,0,x,y);
    }

    inline double bb_f4_dxy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bxy(order,i,j,0,x,y);
    }

    // number 5
    inline double bb_f5(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return B(order,0,j,k,x,y);
    }

    inline double bb_f5_dx(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bx(order,0,j,k,x,y);
    }

    inline double bb_f5_dy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return By(order,0,j,k,x,y);
    }

    inline double bb_f5_dxx(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bxx(order,0,j,k,x,y);
    }

    inline double bb_f5_dyy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Byy(order,0,j,k,x,y);
    }

    inline double bb_f5_dxy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bxy(order,0,j,k,x,y);
    }

    // number 6
    inline double bb_f6(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return B(order,i,0,k,x,y);
    }

    inline double bb_f6_dx(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bx(order,i,0,k,x,y);
    }

    inline double bb_f6_dy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return By(order,i,0,k,x,y);
    }

    inline double bb_f6_dxx(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bxx(order,i,0,k,x,y);
    }

    inline double bb_f6_dyy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Byy(order,i,0,k,x,y);
    }

    inline double bb_f6_dxy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bxy(order,i,0,k,x,y);
    }


  // ORDER 2 & higher
    // Edge functions

    // number 7
    inline double bb_f7(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return B(order,j,i,0,x,y);
    }

    inline double bb_f7_dx(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bx(order,j,i,0,x,y);
    }

    inline double bb_f7_dy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return By(order,j,i,0,x,y);
    }

    inline double bb_f7_dxx(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bxx(order,j,i,0,x,y);
    }

    inline double bb_f7_dyy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Byy(order,j,i,0,x,y);
    }

    inline double bb_f7_dxy(double x, double y, int order=0)
    {
					int i = order-1; 
					int j= 0;
					if(i>0) j = order-i;					
				 return Bxy(order,j,i,0,x,y);
    }

    // number 8
    inline double bb_f8(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return B(order,0,k,j,x,y);
    }

    inline double bb_f8_dx(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bx(order,0,k,j,x,y);
    }

    inline double bb_f8_dy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return By(order,0,k,j,x,y);
    }

    inline double bb_f8_dxx(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bxx(order,0,k,j,x,y);
    }

    inline double bb_f8_dyy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Byy(order,0,k,j,x,y);
    }

    inline double bb_f8_dxy(double x, double y, int order=0)
    {
					int j = order-1; 
					int k =0;
					if(j>0) k = order-j;					
				 return Bxy(order,0,k,j,x,y);
    }

    // number 9
    inline double bb_f9(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return B(order,k,0,i,x,y);
    }

    inline double bb_f9_dx(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bx(order,k,0,i,x,y);
    }

    inline double bb_f9_dy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return By(order,k,0,i,x,y);
    }

    inline double bb_f9_dxx(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bxx(order,k,0,i,x,y);
    }

    inline double bb_f9_dyy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Byy(order,k,0,i,x,y);
    }

    inline double bb_f9_dxy(double x, double y, int order=0)
    {
					int k = order-1; 
					int i =0;
					if(k>0) i = order-k;					
				 return Bxy(order,k,0,i,x,y);
    }


  // ORDER 3 & higher
    // bubble functions, order 3

    // number 10
    inline double bb_f10(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return B(order,i,j,k,x,y);
    }

    inline double bb_f10_dx(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return Bx(order,i,j,k,x,y);
    }

    inline double bb_f10_dy(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return By(order,i,j,k,x,y);
    }

    inline double bb_f10_dxx(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return Bxx(order,i,j,k,x,y);
    }

    inline double bb_f10_dyy(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return Byy(order,i,j,k,x,y);
    }

    inline double bb_f10_dxy(double x, double y, int order=0)
    {
					int i = order - 2; 
					int j =0; 
					if(i>0) j = order-2;
					int k = order-i-j;				
				 return Bxy(order,i,j,k,x,y);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////

    static ShapesetBB::shape_fn_bb bb_tri_fn[] =
    {
      bb_f1,    bb_f2,    bb_f3, 
		   bb_f4,    bb_f5,    bb_f6,
				bb_f7, 		bb_f8, 		bb_f9, 	bb_f10
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dx[] =
    {
      bb_f1_dx,   bb_f2_dx,   bb_f3_dx,  
			 bb_f4_dx,   bb_f5_dx,   bb_f6_dx,
				bb_f7_dx, 		bb_f8_dx, 		bb_f9_dx, 	bb_f10_dx
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dy[] =
    {
      bb_f1_dy,   bb_f2_dy,   bb_f3_dy,   
				bb_f4_dy,   bb_f5_dy,   bb_f6_dy,
					bb_f7_dy, 		bb_f8_dy, 		bb_f9_dy, 	bb_f10_dy
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dxx[] =
    {
      bb_f1_dxx,   bb_f2_dxx,   bb_f3_dxx,
		   bb_f4_dxx,   bb_f5_dxx,   bb_f6_dxx,
				bb_f7_dxx, 		bb_f8_dxx, 		bb_f9_dxx, 	bb_f10_dxx
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dyy[] =
    {
      bb_f1_dyy,   bb_f2_dyy,   bb_f3_dyy,   
				bb_f4_dyy,   bb_f5_dyy,   bb_f6_dyy,
					bb_f7_dyy, 		bb_f8_dyy, 		bb_f9_dyy, 	bb_f10_dyy
    };

    static ShapesetBB::shape_fn_bb bb_tri_fn_dxy[] =
    {
      bb_f1_dxy,   bb_f2_dxy,   bb_f3_dxy,   
				bb_f4_dxy,   bb_f5_dxy,   bb_f6_dxy,
					bb_f7_dxy, 		bb_f8_dxy, 		bb_f9_dxy, 	bb_f10_dxy
    };
////////////
    static int bb_tri_bubble_indices_all_orders[] =
    {
      9

    };
    static int* bb_tri_bubble_indices[4] =
    {
 				NULL, NULL, NULL,
      bb_tri_bubble_indices_all_orders
    };

    static int bb_tri_bubble_count[4] = { 0, 0, 0, 1};

    static int bb_tri_edge_indices_0[8] =  { 0, 1, 1, 0, 3, 6,	6, 3 	};
    static int bb_tri_edge_indices_1[8] =  { 1, 2, 2, 1, 4, 7,	7, 4	};
    static int bb_tri_edge_indices_2[8] =  { 2, 0, 0, 2, 5, 8,	8, 5 };

    static int* bb_tri_edge_indices[3] =
    {
      bb_tri_edge_indices_0,
      bb_tri_edge_indices_1,
      bb_tri_edge_indices_2
    };

    static int bb_tri_vertex_indices[3] = { 0, 1, 2 };

  /*  static int bb_tri_index_to_order[10] =
    {
      1, 1, 1, 
			2, 2, 2, 
			2, 2, 2,	3  //wg. symmetrie 6,7,8 polygrad 2
    };*/

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

 /*   static int* bb_index_to_order[2] =
    {
      bb_tri_index_to_order,
      NULL
    };*/

    int ShapesetBB::get_max_index(ElementMode2D mode) { 
				if(mode == HERMES_MODE_TRIANGLE) return (no_shape_fn-1);
				else return 0.;
			}

    ShapesetBB::ShapesetBB(int order)
    {
			this->space_order = order;
 			no_shape_fn = faculty(2+order)/faculty(order);

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
      //index_to_order = bb_index_to_order;

			index_to_order = new int*[2];
			index_to_order[0] = new int [10];
			index_to_order[1] = NULL;
			for(int i = 0; i< no_shape_fn; i++)
				index_to_order[0][i] = order;


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

      max_order = 3;
      min_order = 1;
      num_components = 1;

      ebias = 2;

      comb_table = NULL;
    }

ShapesetBB::~ShapesetBB(){
		delete index_to_order[0];
		delete index_to_order;
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



    //const int ShapesetBB::max_index[2] = { 9, 0 };
  }
}
