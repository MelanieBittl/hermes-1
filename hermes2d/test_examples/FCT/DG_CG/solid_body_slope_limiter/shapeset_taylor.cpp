#include "shapeset_taylor.h" 


 void AuxiliarySolution::derivatives(double x, double y, double& dx, double& dy) const
  {     
		dx=1.; dy=1.;	
  };

 double AuxiliarySolution::value(double x, double y) const 
 {
	if(x_component)
       return (x-c);
   else 
   	return (y-c);
};

 Ord AuxiliarySolution::ord(double x, double y)  const 
 {
      return Ord(10);
};

 MeshFunction<double>* AuxiliarySolution::clone() const
    {
			return new AuxiliarySolution(this->mesh, x_component, c);

    }


	static double taylor_1(double x, double y, double x_c=NULL, double y_c=NULL,  double rest=NULL)
	{
		return 1.0;
	}
	static double taylor_0(double x, double y, double x_c=NULL, double y_c=NULL,  double rest=NULL)
	{
		return 0.;
	}
	static double taylor_2(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		double a = (x-x_c);
		return a;
	}
	static double taylor_2_dx(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		return 1.;
	}	
	static double taylor_3(double x, double y, double x_c, double y_c,double rest =NULL)
	{
		return (y-y_c);
	}
	static double taylor_3_dy(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		return 1.;
	}	
	static double taylor_4(double x, double y, double x_c, double y_c, double rest)
	{
		double a = (x-x_c)*(x-x_c);
		double b =	rest; //vol_av->calculate_Volume_Average(active_element,0, sln_x, sln_x);
		double d= 2.;	
		return a/d - b/d;
	}
	static double taylor_4_dx(double x, double y, double x_c, double y_c,  double rest=NULL)
	{
		double a = (x-x_c);	
		return a;
	}	
	static double taylor_4_dxx(double x, double y, double x_c, double y_c, double rest=NULL)
	{		
		return 1. ;
	}
	
	static double taylor_5(double x, double y, double x_c, double y_c, double rest)
	{
		double a = (y-y_c)*(y-y_c);
		double b =	rest; //vol_av->calculate_Volume_Average(active_element,0, sln_y, sln_y);
		double d= 2.;	
		return a/d - b/d;
	}
	static double taylor_5_dy(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		double a = (y-y_c);	
		return a ;
	}
	static double taylor_5_dyy(double x, double y, double x_c, double y_c, double rest=NULL)
	{		
		return 1. ;
	}
	static double taylor_6(double x, double y, double x_c, double y_c,  double rest)
	{
		double a = (x-x_c)*(y-y_c);
		double b =	rest; //vol_av->calculate_Volume_Average(active_element,0, sln_x, sln_y);		
		return a - b;
	}
	
	static double taylor_6_dx(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		double a = (y-y_c);				
		return a;
	}

	static double taylor_6_dy(double x, double y, double x_c, double y_c, double rest=NULL)
	{
		double a = (x-x_c);		
		return a;
	}
	static double taylor_6_dxy(double x, double y, double x_c, double y_c, double rest=NULL)
	{		
		return 1.;
	}

 static int vertex_indices[4] = { -1, -1, -1, -1 };

    static int* taylor_vertex_indices[2] =
    {
     vertex_indices, vertex_indices
    };
    
    static int taylor_quad_edge_indices_0[6] =  { -1, -1, -1, -1, -1, -1  };
    static int taylor_quad_edge_indices_1[6] =  { -1, -1, -1, -1, -1, -1  };
    static int taylor_quad_edge_indices_2[6] =  { -1, -1, -1, -1, -1, -1 };
    static int taylor_quad_edge_indices_3[6] =  { -1, -1, -1, -1, -1, -1  };

    int* taylor_quad_edge_indices[4] =
    {
      taylor_quad_edge_indices_0,
      taylor_quad_edge_indices_1,
      taylor_quad_edge_indices_2,
      taylor_quad_edge_indices_3
    };     

    static int** taylor_edge_indices[2] =
    {
      taylor_quad_edge_indices,
      taylor_quad_edge_indices
    };
    

    
    
    static int qb_0_0[] = { 0, };
    static int qb_0_1[] = { 0, 1, };
    static int qb_0_2[] = { 0, 1, 3 };   
    static int qb_1_0[] = { 0, 2, };
    static int qb_1_1[] = { 0, 1, 2 };
    static int qb_1_2[] = { 0, 1, 2, 3}; 
    static int qb_2_0[] = { 0, 2, 4 };
    static int qb_2_1[] = { 0, 1, 2, 4 };
    static int qb_2_2[] = { 0, 1, 2, 3, 4, 5};

    #define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL
	 #define NULL08 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
    int* taylor_quad_bubble_indices[] =
    {
      qb_0_0,  qb_0_1,  qb_0_2,  NULL08,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_1_0,  qb_1_1,  qb_1_2,  NULL08,   NULL, NULL, NULL, NULL, NULL, NULL16,
      qb_2_0,  qb_2_1,  qb_2_2,  NULL08,   NULL, NULL, NULL, NULL, NULL, NULL16,
      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,
      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16,      NULL16, NULL16

    };
    
    static int** taylor_bubble_indices[2] =
    {
      NULL,
      taylor_quad_bubble_indices
    };

    #define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    int taylor_quad_bubble_count[] =
    {
      1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
      2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
      3,  4,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
		zero16,zero16,		zero16,zero16,		zero16,zero16,		zero16,zero16,
		zero16,zero16,		zero16,zero16,		zero16,zero16,		zero16,zero16
    };
    

    static int* taylor_bubble_count[2] =
    {
    	 NULL,
      taylor_quad_bubble_count
    };



    #define oo H2D_MAKE_QUAD_ORDER
    #define XX(a, b) oo(a, b), oo(a, b)

    int taylor_quad_index_to_order[] =
    {
      oo(0, 0),   oo(0, 1),   oo(1, 0),   oo(2, 0),   oo(0, 2),   oo(2, 2)  
     
    };

    int taylor_tri_index_to_order[] =
    {
      0, 
      1, 1,
      2, 2, 2
    };

     int* taylor_index_to_order[2] =
    {
      taylor_tri_index_to_order,
      taylor_quad_index_to_order
    };

		static TaylorShapeset::shape_fn_t taylor_quad_fn[] = 
			{taylor_1, taylor_2, taylor_3, taylor_4, taylor_5, taylor_6};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table[1] = { taylor_quad_fn };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table[2] = { NULL, taylor_quad_shape_fn_table  };

		static TaylorShapeset::shape_fn_t taylor_quad_fn_dx[] =
			 {taylor_0, taylor_2_dx, taylor_0, taylor_4_dx, taylor_0, taylor_6_dx};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table_dx[1] = { taylor_quad_fn_dx };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table_dx[2] = { NULL, taylor_quad_shape_fn_table_dx  };
		
		static TaylorShapeset::shape_fn_t taylor_quad_fn_dy[] = 
				{taylor_0, taylor_0, taylor_3_dy, taylor_0, taylor_5_dy, taylor_6_dy};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table_dy[1] = { taylor_quad_fn_dy };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table_dy[2] = { NULL, taylor_quad_shape_fn_table_dy };

		static TaylorShapeset::shape_fn_t taylor_quad_fn_dxx[] =
			 {taylor_0, taylor_0, taylor_0, taylor_4_dxx, taylor_0, taylor_0};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table_dxx[1] = { taylor_quad_fn_dxx };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table_dxx[2] = { NULL, taylor_quad_shape_fn_table_dxx  };

		static TaylorShapeset::shape_fn_t taylor_quad_fn_dyy[] = 
				{taylor_0, taylor_0, taylor_0, taylor_0, taylor_5_dyy, taylor_0};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table_dyy[1] = { taylor_quad_fn_dyy };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table_dyy[2] = { NULL, taylor_quad_shape_fn_table_dyy };
		
				static TaylorShapeset::shape_fn_t taylor_quad_fn_dxy[] =
			 {taylor_0, taylor_0, taylor_0, taylor_0, taylor_0, taylor_6_dxy};
		TaylorShapeset::shape_fn_t* taylor_quad_shape_fn_table_dxy[1] = { taylor_quad_fn_dxy };
		static TaylorShapeset::shape_fn_t** taylor_shape_fn_table_dxy[2] = { NULL, taylor_quad_shape_fn_table_dxy  };

TaylorShapeset::TaylorShapeset()
{
	active_element =NULL;
//	mesh = NULL;
	//space = NULL;
	vol_av = NULL;
	sln_x = NULL;
	sln_y =NULL;
	
	      max_order = 2;
      num_components = 1;      
      ebias = 2;
      comb_table = NULL;
      
}
TaylorShapeset::TaylorShapeset(Element* e,SpaceSharedPtr<double>  space)
{	
	vol_av = NULL;
	sln_x = NULL;
	sln_y =NULL;
	      max_order = 2;
      num_components = 1;      
      ebias = 2;
      comb_table = NULL;
	init_data(e,space);
}


 void TaylorShapeset::init_data(Element* e,SpaceSharedPtr<double>  new_space)
 {
 			
 			active_element = e;	
 			space = new_space;
 			mesh = space->get_mesh();
 			
 	shape_table[0] = taylor_shape_fn_table;
     shape_table[1] = taylor_shape_fn_table_dx;
     shape_table[2] = taylor_shape_fn_table_dy;
      shape_table[3]= taylor_shape_fn_table_dxx;
      shape_table[4] = taylor_shape_fn_table_dyy;
      shape_table[5] = taylor_shape_fn_table_dxy;


      vertex_indices = taylor_vertex_indices;
      edge_indices = taylor_edge_indices;
      bubble_indices = taylor_bubble_indices;
      bubble_count = taylor_bubble_count;
      index_to_order = taylor_index_to_order;	

      
		int mode = 1;
		x_c =0; y_c = 0;
		
		x_max =e->vn[0]->x; y_max = e->vn[0]->y; 
		x_min =e->vn[0]->x; y_min = e->vn[0]->y;
		if(active_element->is_triangle()) mode = 0;
		for(int i =0; i< active_element->get_nvert();i++)	
		{
			ref_vert[mode][i][0] = e->vn[i]->x;
			ref_vert[mode][i][1] = e->vn[i]->y;
			x_c += active_element->vn[i]->x; 
			y_c += active_element->vn[i]->y;
			if(active_element->vn[i]->x > x_max) x_max = active_element->vn[i]->x;
			else if(active_element->vn[i]->x < x_min) x_min = active_element->vn[i]->x;
			if(active_element->vn[i]->y > y_max) y_max = active_element->vn[i]->y;
			else if(active_element->vn[i]->y < y_min) y_min = active_element->vn[i]->y;
			
		}

 	/*  ref_vert[1][0][0] = -1.0;
      ref_vert[1][0][1] = -1.0;
      ref_vert[1][1][0] =  1.0;
      ref_vert[1][1][1] = -1.0;
      ref_vert[1][2][0] =  1.0;
      ref_vert[1][2][1] =  1.0;
      ref_vert[1][3][0] = -1.0;
      ref_vert[1][3][1] =  1.0;*/


		x_c/= active_element->get_nvert(); 
		y_c /= active_element->get_nvert();

//x_c =0; y_c = 0;
	
	if(vol_av!=NULL) delete vol_av; 
	 if(sln_x!=NULL) delete sln_x; 
	 if(sln_y!=NULL)  delete sln_y;
	 vol_av = NULL;
	sln_x = NULL;
	sln_y =NULL;
	
	   if(vol_av==NULL) vol_av = new VolumeAverage;  	
	   if(sln_x==NULL) sln_x = new AuxiliarySolution(mesh, true, x_c) ;
   	if(sln_y==NULL) sln_y = new AuxiliarySolution(mesh, false, y_c) ;
   	
  rest_taylor_4 = vol_av->calculate_Volume_Average(active_element,0, sln_x, space, sln_x);
   rest_taylor_5 = vol_av->calculate_Volume_Average(active_element,0, sln_y, space, sln_y);
   rest_taylor_6 = vol_av->calculate_Volume_Average(active_element,0, sln_x, space, sln_y);
   
  // rest_taylor_4 = 1./3.;rest_taylor_5 = 1./3.;rest_taylor_6 = 0;  //fuer x_c, y_c  =0 und einheitsquadrat
   
   	delta_x = (x_max - x_min)/2.;
		delta_y = (y_max - y_min)/2.;    

			//x_c =0; y_c = 0;
 }

TaylorShapeset::~TaylorShapeset()
{
	 if(vol_av!=NULL) delete vol_av; 
	 if(sln_x!=NULL) delete sln_x; 
	 if(sln_y!=NULL)  delete sln_y;
	 vol_av = NULL;
	sln_x = NULL;
	sln_y =NULL;

}

double TaylorShapeset::get_value(int n, int index, double x, double y, int component, ElementMode2D mode)
{     
	 if(index >= 0)
      {
        TaylorShapeset::shape_fn_t** shape_expansion = shape_table[n][mode];
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
        	double rest;
        		switch(index)
			    {
			    	case 0: 	rest = 0; break;
			    	case 1:	rest = 0; break;
			    	case 2:  rest = 0; break;
			    	case 3:	rest = rest_taylor_4;break;
			    	case 4:	rest = rest_taylor_5; break;
			    	case 5:	rest = rest_taylor_6;  break;        	
			    	default: 
			    		throw Exceptions::Exception("TaylorSolution::get_fn_values, false idx =, %d",index);			   
			    }	       
      
         return shape_expansion[component][index](x, y, x_c, y_c, rest);
       
        }
      }
      else
        return get_constrained_value(n, index, x, y, component, mode);




}



int TaylorShapeset::get_max_index(ElementMode2D mode) { return max_index[mode]; }
const int TaylorShapeset::max_index[2] = { 5,5 };
