 #include "solution_slopelimiter.h"

//Detector for Oszillations

KuzminOscillationDetector::KuzminOscillationDetector(SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double> solution) : space(space), solution(solution)
{
  mesh = space->get_mesh();
	vol_average = new VolumeAverage;

  int max_id = mesh->get_max_node_id();
	u_min_i = new double[max_id];
	u_max_i = new double[max_id];
	u_dx_min_i = new double[max_id];
	u_dx_max_i = new double[max_id];
	u_dy_min_i = new double[max_id];
	u_dy_max_i = new double[max_id];   

  int max_elem = mesh->get_max_element_id();   
  u_c = new double[max_elem]; 
  u_dx_K = new double[max_elem];  
  u_dy_K= new double[max_elem];  
  u_dxx_K= new double[max_elem]; 
  u_dyy_K= new double[max_elem];   
  u_dxy_K= new double[max_elem];  
 
  alpha_i_first_order = new double[max_elem]; 
  alpha_i_second_order= new double[max_elem];   
   
};

KuzminOscillationDetector::~KuzminOscillationDetector()
{
	
  delete [] u_min_i;
  delete [] u_max_i;
  delete [] u_dx_min_i;
  delete [] u_dx_max_i;
  delete [] u_dy_min_i;
  delete [] u_dy_max_i;  
  delete [] u_c;  
 	delete [] u_dx_K; 
  delete [] u_dy_K;
  delete [] u_dxx_K;
  delete [] u_dyy_K;
  delete [] u_dxy_K;  
  delete [] alpha_i_first_order;
  delete [] alpha_i_second_order;  
  delete vol_average;  

};



  
void KuzminOscillationDetector::get_delta(Element* e, double &delta_x, double &delta_y)
{
	double x_max, x_min, y_max, y_min;
	x_max = e->vn[0]->x; 
	y_max = e->vn[0]->y;
	x_min = e->vn[0]->x; 
	y_min = e->vn[0]->y;
	for(unsigned int i = 1; i < e->get_nvert(); i++)
	{       
		if(e->vn[i]->x > x_max) x_max = e->vn[i]->x;
		if(e->vn[i]->y > y_max) y_max = e->vn[i]->y;
		 if(e->vn[i]->x < x_min) x_min = e->vn[i]->x;
		if(e->vn[i]->y < y_min) y_min = e->vn[i]->y;
	}
	delta_x = (x_max - x_min)/2.; 
	delta_y = (y_max - y_min)/2.;  
	if(delta_x==0) Hermes::Exceptions::Exception("delta_x==0. ");
	if(delta_y==0) Hermes::Exceptions::Exception("delta_y==0. ");  
}



void KuzminOscillationDetector::compute_all(bool only_p1)
{
		Element* e; Node* vn;
	double delta_x, delta_y, x_max, x_min, y_max, y_min;
   for_all_vertex_nodes(vn, mesh)
   {   
   	u_min_i[vn->id] =	std::numeric_limits<double>::infinity();
   	u_max_i[vn->id] = -std::numeric_limits<double>::infinity();
   	u_dx_min_i[vn->id] =	std::numeric_limits<double>::infinity();
   	u_dx_max_i[vn->id] = -std::numeric_limits<double>::infinity();	
   	u_dy_min_i[vn->id] =	std::numeric_limits<double>::infinity();
   	u_dy_max_i[vn->id] = -std::numeric_limits<double>::infinity();		
   }
		
		bool set_delta= false;
    for_all_active_elements(e, this->mesh)
    {
			alpha_i_first_order[e->id] =1.;
		   alpha_i_second_order[e->id]=1.;		
			if(set_delta==false)
			{   
				get_delta(e, delta_x, delta_y);
			  set_delta = true;
		  }     

		   
		   double u=0.; double u_1=0.; double u_2=0.;double u_3=0.;
			find_centroid_values(e, u);
			 u_c[e->id] = u;
			 find_centroid_derivatives(e, u_1, u_2);
			 u_dx_K[e->id] = delta_x*u_1;
			 u_dy_K[e->id] = delta_y*u_2;
			 find_second_centroid_derivatives(e, u_1, u_2,u_3);  	 

		 u_dxx_K[e->id] = (delta_x*delta_x)*u_1;
		 u_dxy_K[e->id] = (delta_x*delta_y)*u_2;
		 u_dyy_K[e->id] = (delta_y*delta_y)*u_3;
			
			for(int i =0; i<e->get_nvert(); i++)
			{
				vn = e->vn[i];			
				if(u_min_i[vn->id] > u_c[e->id]) 
					u_min_i[vn->id] =	u_c[e->id];
				if(u_max_i[vn->id] < u_c[e->id]) 	
					u_max_i[vn->id] = u_c[e->id];
				if(u_dx_min_i[vn->id] > u_dx_K[e->id]) 
					u_dx_min_i[vn->id] =	u_dx_K[e->id];
				if(u_dx_max_i[vn->id] < u_dx_K[e->id]) 	
					u_dx_max_i[vn->id] = u_dx_K[e->id];
				if(u_dy_min_i[vn->id] > u_dy_K[e->id]) 
					u_dy_min_i[vn->id] =	u_dy_K[e->id];
				if(u_dy_max_i[vn->id] < u_dy_K[e->id]) 	
					u_dy_max_i[vn->id] = u_dy_K[e->id];			
			}
			
		}
		
		double u_i;
		double u_dx_i;
		double u_dy_i;
	 	
    for_all_active_elements(e, this->mesh)
      {
			double c_ref_x =0; double c_ref_y=0; 
			 x_max = e->vn[0]->x; 
			 y_max = e->vn[0]->y;
			 x_min = e->vn[0]->x; 
			 y_min = e->vn[0]->y;
      	for(unsigned int i = 0; i < e->get_nvert(); i++)
         {
         	c_ref_x +=e->vn[i]->x; 
         	c_ref_y +=e->vn[i]->y; 

         }
         c_ref_y /=e->get_nvert(); c_ref_x/=e->get_nvert(); 
         
			if(set_delta==false)
			{   
					get_delta(e, delta_x, delta_y);
				   set_delta = true;
		   }   
          
		  	for(unsigned int i = 0; i < e->get_nvert(); i++)
			 {
				vn = e->vn[i];					
     
      		u_i = 	 u_dx_K[e->id]  *(vn->x- c_ref_x)/delta_x + u_dy_K[e->id]  * (vn->y -c_ref_y)/delta_y;
      		u_dx_i =  u_dxx_K[e->id] *(vn->x- c_ref_x)/delta_x + u_dxy_K[e->id] * (vn->y -c_ref_y)/delta_y;
      		u_dy_i =  u_dxy_K[e->id] *(vn->x- c_ref_x)/delta_x	+ u_dyy_K[e->id] * (vn->y -c_ref_y)/delta_y;  				
				
			
				if((u_min_i[vn->id] - u_c[e->id]) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u_i_min[vertex_i] - u_c)>0. ");	
				if((u_max_i[vn->id] - u_c[e->id]) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_i_max[vertex_i] - u_c)<0"); 
		
				if(u_i < 0)
				{
				  if(((u_min_i[vn->id] - u_c[e->id]) / u_i)  < alpha_i_first_order[e->id])
				    alpha_i_first_order[e->id] = (u_min_i[vn->id] - u_c[e->id]) / u_i ;
				}
				else if(u_i > 0)
				{
				  if(((u_max_i[vn->id] - u_c[e->id]) / u_i) < alpha_i_first_order[e->id])
				    alpha_i_first_order[e->id] = (u_max_i[vn->id] - u_c[e->id]) / u_i ;
				}
				
				if(!only_p1)
				{
					if((u_dx_min_i[vn->id] - u_dx_K[e->id]) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u__d_i_min[vertex_i][0] - u__dx_c)>0. ");	
					if((u_dx_max_i[vn->id] - u_dx_K[e->id]) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_di_max[vertex_i][0] - u_dx_c)<0"); 
					if((u_dy_min_i[vn->id] - u_dy_K[e->id]) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u__d_i_min[vertex_i][1] - u__dy_c)>0. ");	
					if((u_dy_max_i[vn->id] - u_dy_K[e->id]) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_di_max[vertex_i][1] - u_dy_c)<0"); 

					// dx.
					if(u_dx_i < 0)
					{
					  if(((u_dx_min_i[vn->id] - u_dx_K[e->id]) / u_dx_i)  < alpha_i_second_order[e->id])
						 alpha_i_second_order[e->id] = (u_dx_min_i[vn->id] - u_dx_K[e->id]) / u_dx_i ;
					}
					else if(u_dx_i > 0)
					{
					  if(((u_dx_max_i[vn->id] - u_dx_K[e->id]) / u_dx_i) < alpha_i_second_order[e->id])
						 alpha_i_second_order[e->id] = (u_dx_max_i[vn->id] - u_dx_K[e->id]) / u_dx_i;
					}
					// dy.
					if(u_dy_i < 0)
					{
					  if(((u_dy_min_i[vn->id] - u_dy_K[e->id]) / u_dy_i) < alpha_i_second_order[e->id])
						 alpha_i_second_order[e->id] = (u_dy_min_i[vn->id] - u_dy_K[e->id]) / u_dy_i;
					}
					else if(u_dy_i > 0)
					{
					  if(((u_dy_max_i[vn->id] - u_dy_K[e->id]) / u_dy_i) < alpha_i_second_order[e->id])
						 alpha_i_second_order[e->id] = (u_dy_max_i[vn->id] - u_dy_K[e->id]) / u_dy_i;
					} 
				}

				
			 }
			 if(!only_p1)
 		 			if(alpha_i_first_order[e->id] < alpha_i_second_order[e->id])
	   				alpha_i_first_order[e->id] = alpha_i_second_order[e->id];
	   	 else
	   	 	alpha_i_second_order[e->id] =0;
	   	
	   	if((alpha_i_first_order[e->id]<0.)||(alpha_i_first_order[e->id]>1.))   Hermes::Exceptions::Exception("alpha_first_order. ");
	   	if((alpha_i_second_order[e->id]<0.)||(alpha_i_second_order[e->id]>1.)) Hermes::Exceptions::Exception("alpha_sec_order. ");
      
      }
		

}




void KuzminOscillationDetector::find_centroid_values(Element* e, double &u_c)
{
    u_c = vol_average->calculate_Volume_Average(e,0, solution, space);

   /* solution->set_active_element(e);  
        if(e->get_mode() == HERMES_MODE_TRIANGLE)
				{
          u_c = (dynamic_cast<Solution<double>*>(solution.get()))->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 0);
				
        }else{
          u_c = (dynamic_cast<Solution<double>*>(solution.get()))->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 0);
				}*/

}

void KuzminOscillationDetector::find_centroid_derivatives(Element* e, double &u_dx_c, double &u_dy_c)
{

    solution->set_active_element(e);  
       if(e->get_mode() == HERMES_MODE_TRIANGLE)
				{
          u_dx_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 1);
					u_dy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 2);
        }else{
          u_dx_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 1);
					u_dy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 2);
				} 

  
}

void KuzminOscillationDetector::find_second_centroid_derivatives(Element* e, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c)
{


    solution->set_active_element(e);  
       if(e->get_mode() == HERMES_MODE_TRIANGLE)
				{
          u_dxx_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 3);
					u_dyy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 4);
					u_dxy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, 5);
        }else{
          u_dxx_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 3);
					u_dyy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 4);
					u_dxy_c = solution.get_solution()->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, 5);
				}   
  
}




    static double3* cheb_tab_tri[11];
    static double3* cheb_tab_quad[11];
    static int      cheb_np_tri[11];
    static int      cheb_np_quad[11];

    static double3** cheb_tab[2] = { cheb_tab_tri, cheb_tab_quad };
    static int*      cheb_np[2]  = { cheb_np_tri,  cheb_np_quad  };

    static class Quad2DCheb : public Quad2D
    {
    public:

      Quad2DCheb()
      {
        max_order[0]  = max_order[1]  = 10;
        num_tables[0] = num_tables[1] = 11;
        tables = cheb_tab;
        np = cheb_np;

        tables[0][0] = tables[1][0] = NULL;
        np[0][0] = np[1][0] = 0;

        int i, j, k, n, m;
        double3* pt;
        for (int mode_i = 0; mode_i <= 1; mode_i++)
        {
          for (k = 0; k <= 10; k++)
          {
            np[mode_i][k] = n = mode_i ? Hermes::sqr(k + 1) : (k + 1)*(k + 2)/2;
            tables[mode_i][k] = pt = new double3[n];

            for (i = k, m = 0; i >= 0; i--)
              for (j = k; j >= (mode_i ? 0 : k-i); j--, m++) {
                pt[m][0] = k ? Hermes::cos(j * M_PI / k) : 1.0;
                pt[m][1] = k ? Hermes::cos(i * M_PI / k) : 1.0;
                pt[m][2] = 1.0;
              }
          }
        }
      };

      ~Quad2DCheb()
      {
        for (int mode_i = 0; mode_i <= 1; mode_i++)
          for (int k = 1; k <= 10; k++)
            delete [] tables[mode_i][k];
      }

      virtual void dummy_fn() {}
    } g_quad_2d_cheb;
    
    
    
        static struct mono_lu_init
    {
    public:

      // this is a set of LU-decomposed matrices shared by all Solutions
      double** mat[2][11];
      int* perm[2][11];

      mono_lu_init()
      {
        memset(mat, 0, sizeof(mat));
      }

      ~mono_lu_init()
      {
        for (int m = 0; m <= 1; m++)
          for (int i = 0; i <= 10; i++)
            if(mat[m][i] != NULL) {
              delete [] mat[m][i];
              delete [] perm[m][i];
            }
      }
    }
    mono_lu;
  
  
  
  
  
  SlopeLimiterSolution::SlopeLimiterSolution(MeshFunctionSharedPtr<double> ref_sln, SpaceSharedPtr<double> space): Solution<double>(space->get_mesh()), space(space), ref_sln(ref_sln) 
   {  	
			detector = new KuzminOscillationDetector(space, ref_sln); 			
  };

  
   
 void SlopeLimiterSolution::limit_solution_according_to_detector(bool p1_only)
    {
     int o;			
      free();
      this->space_type = space->get_type();

		if(this->space_type != HERMES_L2_SPACE)
			throw Exceptions::Exception("TaylorSolution need to be in L2 Space.");
      this->num_components = ref_sln->get_num_components();
      sln_type = HERMES_SLN;		
      
      // Copy the mesh.
      this->mesh = space->get_mesh();

      // Allocate the coefficient arrays.
      num_elems = this->mesh->get_max_element_id();
      if(elem_orders != NULL)
        delete [] elem_orders;
      elem_orders = new int[num_elems];
      memset(elem_orders, 0, sizeof(int) * num_elems);
      for (int l = 0; l < this->num_components; l++)
      {
        if(elem_coeffs[l] != NULL)
          delete [] elem_coeffs[l];
        elem_coeffs[l] = new int[num_elems];
        memset(elem_coeffs[l], 0, sizeof(int) * num_elems);
      }

      // Obtain element orders, allocate mono_coeffs.
      Element* e;
      num_coeffs = 0;
      for_all_active_elements(e, this->mesh)
      {
        this->mode = e->get_mode();
        o = space->get_element_order(e->id);
        o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
        for (unsigned int k = 0; k < e->get_nvert(); k++)
        {
          int eo = space->get_edge_order(e, k);
          if(eo > o) o = eo;
        }

        num_coeffs += this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;
        elem_orders[e->id] = o;
      }
      num_coeffs *= this->num_components;
      if(mono_coeffs != NULL)
        delete [] mono_coeffs;
      mono_coeffs = new double[num_coeffs];

      // Express the solution on elements as a linear combination of monomials.

   		 double* mono = mono_coeffs;
      TaylorShapeset *shapeset_taylor = new TaylorShapeset;
	    Shapeset* old_shapeset = space->get_shapeset();
      MeshSharedPtr ref_mesh(new Mesh);
				ref_mesh->copy(this->mesh);
			Space<double>* ref_space = new L2Space<double>(ref_mesh, 2);  

      for_all_active_elements(e, ref_mesh)
      {
        ref_space->set_element_order_internal(e->id, space->get_element_order(e->id));
      }
     
       Quad2D* quad = &g_quad_2d_cheb;
      double* alpha = new double[2];
      double u_c, u_dx_c, u_dy_c, u_dxx_c, u_dxy_c, u_dyy_c;  
      u_c =0;
      u_dy_c =0;
      u_dx_c=0;
      u_dxx_c=0;
      u_dxy_c=0;
      u_dyy_c =0;      
	       
      detector->compute_all(p1_only);      
     refmap->set_quad_2d(quad);
      AsmList<double> al;
      for_all_active_elements(e, this->mesh)
      {
     		 shapeset_taylor->init_data(e, space);
          
	 		 ref_space->set_shapeset(shapeset_taylor);
	     	 ref_space->assign_dofs();
	
	 		alpha[0] = detector->get_alpha_first_order()[e->id];
			 alpha[1] = detector->get_alpha_second_order()[e->id];  
  			
  			detector->find_centroid_values(e, u_c);
  			detector->find_centroid_derivatives(e, u_dx_c, u_dy_c);
 			detector->find_second_centroid_derivatives(e, u_dxx_c, u_dxy_c, u_dyy_c);
 						
			//alpha[0] = 1.; 
			//alpha[1] = 1.;
			
        this->mode = e->get_mode();
        o = elem_orders[e->id];
        int np = quad->get_num_points(o, e->get_mode());
        
        ref_space->get_element_assembly_list(e, &al);            						
						
		  refmap->set_active_element(e);
			double* x =  refmap->get_phys_x(o);
			double* y =  refmap->get_phys_y(o);		

       for (int l = 0; l < this->num_components; l++)
        {
          // Obtain solution values for the current element.
          double* val = mono;
          elem_coeffs[l][e->id] = (int) (mono - mono_coeffs);
          memset(val, 0, sizeof(double)*np);
          for (unsigned int k = 0; k < al.get_cnt(); k++)
          {
          if(al.get_idx()[k]==-1) continue; //Bei l2_semi_cg!!!!
				double u_new=0;
				switch(al.get_idx()[k])
				{
					case 0: 	u_new = u_c;   break;
					case 1:	u_new = u_dx_c*alpha[0]; break;
					case 2:  u_new = u_dy_c*alpha[0]; break;
					case 3:	u_new = u_dxx_c*alpha[1];break;
					case 4:	u_new = u_dyy_c*alpha[1];break;
					case 5:	u_new = u_dxy_c*alpha[1];break;        	
					default: 
					throw Exceptions::Exception("TaylorSolution::set_coeff_vector(), false idx =, %d",al.get_idx()[k]);			   
				}	
			  

				double coef = al.get_coef()[k] * u_new;  
				double3* pts  = quad->get_points(o, e->get_mode());				


				for (int i = 0; i < np; i++)
				{			
					val[i] += pts[i][2]*shapeset_taylor->get_value(0, al.get_idx()[k], x[i], y[i], 0, e->get_mode())* coef;				
				}        
          }
          mono += np;

          // solve for the monomial coefficients
          if(mono_lu.mat[this->mode][o] == NULL)
            mono_lu.mat[this->mode][o] = calc_mono_matrix(o, mono_lu.perm[this->mode][o]);
          lubksb(mono_lu.mat[this->mode][o], np, mono_lu.perm[this->mode][o], val);
        }

      }

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL.\n");
      init_dxdy_buffer();
      this->element = NULL;
       
      	delete ref_space;
        delete [] alpha;
        delete shapeset_taylor;
       
    }
  
  
    
void SlopeLimiterSolution::set_diff_linear()
  {     int o;			
	      free();
      this->space_type = space->get_type();
      
		if(this->space_type != HERMES_L2_SPACE)
			throw Exceptions::Exception("TaylorSolution need to be in L2 Space.");
      this->num_components = ref_sln->get_num_components();
      sln_type = HERMES_SLN;		
      
      // Copy the mesh.
      this->mesh = space->get_mesh();

      // Allocate the coefficient arrays.
      num_elems = this->mesh->get_max_element_id();
      if(elem_orders != NULL)
        delete [] elem_orders;
      elem_orders = new int[num_elems];
      memset(elem_orders, 0, sizeof(int) * num_elems);
      for (int l = 0; l < this->num_components; l++)
      {
        if(elem_coeffs[l] != NULL)
          delete [] elem_coeffs[l];
        elem_coeffs[l] = new int[num_elems];
        memset(elem_coeffs[l], 0, sizeof(int) * num_elems);
      }

      // Obtain element orders, allocate mono_coeffs.
      Element* e;
      num_coeffs = 0;
      for_all_active_elements(e, this->mesh)
      {
        this->mode = e->get_mode();
        o = space->get_element_order(e->id);
        o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
        for (unsigned int k = 0; k < e->get_nvert(); k++)
        {
          int eo = space->get_edge_order(e, k);
          if(eo > o) o = eo;
        }

        num_coeffs += this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;
        elem_orders[e->id] = o;
      }
      num_coeffs *= this->num_components;
      if(mono_coeffs != NULL)
        delete [] mono_coeffs;
      mono_coeffs = new double[num_coeffs];

      // Express the solution on elements as a linear combination of monomials.

    double* mono = mono_coeffs;
      TaylorShapeset *shapeset_taylor = new TaylorShapeset;
	    Shapeset* old_shapeset = space->get_shapeset();
      MeshSharedPtr ref_mesh(new Mesh);
      Space<double>* ref_space =NULL;
      if(old_shapeset->get_space_type()==HERMES_L2_SPACE)
      {
      		ref_space = new L2Space<double>();      		
      }
      else if(old_shapeset->get_space_type() == HERMES_H1_SPACE)
      {
      		 ref_space = new L2_SEMI_CG_Space<double>();
      }		 
      else
      	throw Exceptions::Exception("TaylorSolution:: shapeset isn't in L2 or L2_SEMI_CG.");	
      ref_space->copy(space, ref_mesh);	     
     
      Quad2D* quad = &g_quad_2d_cheb;
      double* alpha = new double[2];
      double u_c, u_dx_c, u_dy_c;  
      u_c =0;
      u_dy_c =0;
      u_dx_c=0;
                 
      detector->compute_all(true);      
     refmap->set_quad_2d(quad);
      
      for_all_active_elements(e, this->mesh)
      {
     		 shapeset_taylor->init_data(e, space);

	 		 ref_space->set_shapeset(shapeset_taylor);
	     	 ref_space->assign_dofs();			
	 		alpha[0] = detector->get_alpha_first_order()[e->id];
	 		  			
  			detector->find_centroid_values(e, u_c);
  			detector->find_centroid_derivatives(e, u_dx_c, u_dy_c);

			
        this->mode = e->get_mode();
        o = elem_orders[e->id];
        int np = quad->get_num_points(o, e->get_mode());

        AsmList<double> al;
        ref_space->get_element_assembly_list(e, &al);
           						
						
		  refmap->set_active_element(e);
			double* x =  refmap->get_phys_x(o);
			double* y =  refmap->get_phys_y(o);		

       for (int l = 0; l < this->num_components; l++)
        {
          // Obtain solution values for the current element.
          double* val = mono;
          elem_coeffs[l][e->id] = (int) (mono - mono_coeffs);
          memset(val, 0, sizeof(double)*np);
          for (unsigned int k = 0; k < al.get_cnt(); k++)
          {
          if(al.get_idx()[k]==-1) continue; //Bei l2_semi_cg!!!!
				double u_new=0;
				switch(al.get_idx()[k])
				{
					case 0: 	u_new = 0;   break;
					case 1:	u_new =  u_dx_c*(1.-alpha[0]); break;
					case 2:  u_new =  u_dy_c*(1.-alpha[0]); break;
					case 3:	u_new = 0;break;
					case 4:	u_new = 0;break;
					case 5:	u_new = 0;break;        	
					default: 
					throw Exceptions::Exception("TaylorSolution::set_coeff_vector(), false idx =, %d",al.get_idx()[k]);			   
				}	
			  
   

				double coef = al.get_coef()[k] * u_new;  
				double3* pts  = quad->get_points(o, e->get_mode());	

				for (int i = 0; i < np; i++)			
						val[i] += pts[i][2]*shapeset_taylor->get_value(0, al.get_idx()[k], x[i], y[i], 0, e->get_mode())* coef;					
 
				        
          }
          mono += np;

          // solve for the monomial coefficients
          if(mono_lu.mat[this->mode][o] == NULL)
            mono_lu.mat[this->mode][o] = calc_mono_matrix(o, mono_lu.perm[this->mode][o]);
          lubksb(mono_lu.mat[this->mode][o], np, mono_lu.perm[this->mode][o], val);
        }
 
      }

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL.\n");
      init_dxdy_buffer();
      this->element = NULL;
       
      
        delete [] alpha;
        delete ref_space;
        delete shapeset_taylor;
       
    }
  

    
  
