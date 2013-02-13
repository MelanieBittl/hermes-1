 #include "solution_slopelimiter.h"

//Detector for Oszillations

KuzminOscillationDetector::KuzminOscillationDetector(Space<double>* space, Solution<double>* solution) : space(space), solution(solution)
{
  mesh = space->get_mesh();
   VolumeAverage* vol_average = new VolumeAverage;
   int max_id=mesh->get_max_node_id();
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


  double KuzminOscillationDetector::get_alpha_first_order(Element* e)
  {
    if(e->get_nvert() == 3)
    	throw Hermes::Exceptions::Exception("So far this limiter is implemented just for quads.");
    double u_c, u_dx_c, u_dy_c;
    find_centroid_values(e, u_c);
    find_centroid_derivatives(e, u_dx_c, u_dy_c);

    // Vertex values.
    double u_i[4];
    find_vertex_values(e, u_i);

    // Boundaries for alpha_i calculation.
    double u_i_min_first_order[4];
    double u_i_max_first_order[4];    
      for(int j = 0; j < 4; j++)
      {
        u_i_min_first_order[j] = u_c;
        u_i_max_first_order[j] = u_c;
      }
    find_u_i_min_max_first_order(e, u_i_min_first_order, u_i_max_first_order);

    // alpha_i calculation.
    double alpha_i_first_order;
    find_alpha_i_first_order(u_i_min_first_order, u_i_max_first_order, u_c, u_i, alpha_i_first_order);   
      if(alpha_i_first_order <0.) Hermes::Exceptions::Exception(" KuzminOscillationDetector::get_alpha_all_ordersalpha_i_first_order <0. ");	
  	return alpha_i_first_order;
  
  }
  double KuzminOscillationDetector::get_alpha_second_order(Element* e)
  {
      if(e->get_nvert() == 3)
      throw Hermes::Exceptions::Exception("So far this limiter is implemented just for quads.");
    double u_dx_c, u_dy_c, u_dxx_c, u_dxy_c, u_dyy_c;
    find_centroid_derivatives(e, u_dx_c, u_dy_c);
    find_second_centroid_derivatives(e, u_dxx_c, u_dxy_c, u_dyy_c);

    // Vertex values.
    double u_d_i[4][2];
    find_vertex_derivatives(e, u_d_i);
    
    // Boundaries for alpha_i calculation.
    double u_d_i_min_second_order[4][2];
    double u_d_i_max_second_order[4][2];

      for(int j = 0; j < 4; j++)
      {
          u_d_i_min_second_order[j][0] = u_dx_c;
          u_d_i_max_second_order[j][0] = u_dx_c;
          u_d_i_min_second_order[j][1] = u_dy_c;
          u_d_i_max_second_order[j][2] = u_dy_c;
       }
    find_u_i_min_max_second_order(e, u_d_i_min_second_order, u_d_i_max_second_order);

    // alpha_i calculation.
    double alpha_i_second_order;
    find_alpha_i_second_order(u_d_i_min_second_order, u_d_i_max_second_order, u_dx_c, u_dy_c, u_d_i, alpha_i_second_order);
      if(alpha_i_second_order<0.) Hermes::Exceptions::Exception(" KuzminOscillationDetector::get_alpha_all_ordersalpha_i_second_order <0. ");	
  	 return alpha_i_second_order;  
  }

void KuzminOscillationDetector::get_alpha_all_orders(Element* e, double* alpha,double &u_c, double &u_dx_c, double &u_dy_c, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c)
{
  if(e->get_nvert() == 3)
      throw Hermes::Exceptions::Exception("So far this limiter is implemented just for quads.");
    //double u_c, u_dx_c, u_dy_c, u_dxx_c, u_dxy_c, u_dyy_c;
    find_centroid_values(e, u_c);
    find_centroid_derivatives(e, u_dx_c, u_dy_c);
    find_second_centroid_derivatives(e, u_dxx_c, u_dxy_c, u_dyy_c);
    
  ////--------------------------  
  double delta_x, delta_y, x_max, x_min, y_max, y_min;
    		double c_ref_x =0; double c_ref_y=0; 
    			x_max = e->vn[0]->x; 
				y_max = e->vn[0]->y;
				x_min = e->vn[0]->x; 
				y_min = e->vn[0]->y;
		   	for(unsigned int i = 0; i < 4; i++)
		      {      c_ref_x +=e->vn[i]->x; c_ref_y +=e->vn[i]->y;  
		      	if(e->vn[i]->x > x_max) x_max = e->vn[i]->x;
		      	if(e->vn[i]->y > y_max) y_max = e->vn[i]->y;
		      	 if(e->vn[i]->x < x_min) x_min = e->vn[i]->x;
		      	if(e->vn[i]->y < y_min) y_min = e->vn[i]->y;
		      }
		      delta_x = (x_max - x_min)/2.; 
		      delta_y = (y_max - y_min)/2.;
		      c_ref_y /=4.; c_ref_x/=4.; 
		 /*   u_dx_c *= delta_x;
		 u_dy_c *= delta_y;
		 u_dxx_c *= (delta_x*delta_x);
		 u_dxy_c *= (delta_x*delta_y);
		 u_dyy_c *= (delta_y*delta_y);*/
	///--------------------------------------	      

    // Vertex values.
    double u_i[4];
   // find_vertex_values(e, u_i);
    double u_d_i[4][2];
   // find_vertex_derivatives(e, u_d_i);
    
    
    //--------------
    Node* vn;
    for(int i=0;i<e->get_nvert();i++)
    {
    		vn=e->vn[i];
          /*  u_i[i] = 		u_c+ u_dx_c  *(vn->x- c_ref_x)/delta_x + u_dy_c  * (vn->y -c_ref_y)/delta_y;
      		u_d_i[i][0] = u_dx_c + u_dxx_c *(vn->x- c_ref_x)/delta_x + u_dxy_c * (vn->y -c_ref_y)/delta_y;
      		u_d_i[i][1] = u_dy_c + u_dxy_c *(vn->x- c_ref_x)/delta_x	+ u_dyy_c * (vn->y -c_ref_y)/delta_y; */
      		
      		u_i[i] = 		u_c+ u_dx_c  *(vn->x- c_ref_x) + u_dy_c  * (vn->y -c_ref_y);
      		u_d_i[i][0] = u_dx_c + u_dxx_c *(vn->x- c_ref_x)+ u_dxy_c * (vn->y -c_ref_y);
      		u_d_i[i][1] = u_dy_c + u_dxy_c *(vn->x- c_ref_x)	+ u_dyy_c * (vn->y -c_ref_y);
     }	
    
    
    
    //----------
    
    

    // Boundaries for alpha_i calculation.
    double u_i_min_first_order[4];
    double u_i_max_first_order[4];
    double u_d_i_min_second_order[4][2];
    double u_d_i_max_second_order[4][2];    
      for(int j = 0; j < 4; j++)
      {
        u_i_min_first_order[j] = u_c;
        u_i_max_first_order[j] = u_c;

          u_d_i_min_second_order[j][0] = u_dx_c;
          u_d_i_max_second_order[j][0] = u_dx_c;
          u_d_i_min_second_order[j][1] = u_dy_c;
          u_d_i_max_second_order[j][1] = u_dy_c;
       
      }
    find_u_i_min_max_first_order(e, u_i_min_first_order, u_i_max_first_order);
    find_u_i_min_max_second_order(e, u_d_i_min_second_order, u_d_i_max_second_order);

    // alpha_i calculation.
    double alpha_i_first_order;
    find_alpha_i_first_order(u_i_min_first_order, u_i_max_first_order, u_c, u_i, alpha_i_first_order);
    // alpha_i calculation.
    double alpha_i_second_order;
    find_alpha_i_second_order(u_d_i_min_second_order, u_d_i_max_second_order, u_dx_c, u_dy_c, u_d_i, alpha_i_second_order); 
    
    if((alpha_i_first_order <0.)||(alpha_i_first_order >1.)) Hermes::Exceptions::Exception(" KuzminOscillationDetector::get_alpha_all_ordersalpha_i_first_order <0. ");	
    if((alpha_i_second_order<0.)||(alpha_i_second_order>1.) ) Hermes::Exceptions::Exception(" KuzminOscillationDetector::get_alpha_all_ordersalpha_i_second_order <0. ");	
    
 	 if(alpha_i_first_order > alpha_i_second_order)
 		alpha[0] = alpha_i_first_order;
 	  else 
 		 alpha[0] = alpha_i_second_order; 	
 	
    	alpha[1] = alpha_i_second_order;	
}


  double KuzminOscillationDetector::get_delta_x(Element* e)
  {
      	double  x_max, x_min;    	
    			x_max = e->vn[0]->x; 		
				x_min = e->vn[0]->x; 
				
		   	for(unsigned int i = 1; i < 4; i++)
		      {     
		      	if(e->vn[i]->x > x_max) x_max = e->vn[i]->x;		      	
		      	 if(e->vn[i]->x < x_min) x_min = e->vn[i]->x;	      	
		      }
		  return (x_max-x_min)/2.;
  
  }
  double KuzminOscillationDetector::get_delta_y(Element* e)
  {
   		double  y_max, y_min;    	
    			y_max = e->vn[0]->y; 		
				y_min = e->vn[0]->y; 				
		   	for(unsigned int i = 1; i < 4; i++)
		      {     
		      	if(e->vn[i]->y > y_max) y_max = e->vn[i]->y;		      	
		      	 if(e->vn[i]->y < y_min) y_min = e->vn[i]->y;	      	
		      }
		  return (y_max-y_min)/2.;
  
  }

void KuzminOscillationDetector::compute_all()
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
   /*	u_min_i[vn->id] =	100000000000000000;
   	u_max_i[vn->id] = -1000000000000000000000;
   	u_dx_min_i[vn->id] =	10000000000000000000000000;
   	u_dx_max_i[vn->id] = -10000000000000000000000000000000;	
   	u_dy_min_i[vn->id] =	100000000000000000000000000000000;
   	u_dy_max_i[vn->id] = -1000000000000000000000000000000000000000000;*/
   }
		
		bool set_delta= false;
      for_all_active_elements(e, this->mesh)
      {
			alpha_i_first_order[e->id] =1.;
		   alpha_i_second_order[e->id]=1.;
		
			if(set_delta==false)
			{   
					x_max = e->vn[0]->x; 
					y_max = e->vn[0]->y;
					x_min = e->vn[0]->x; 
					y_min = e->vn[0]->y;
					for(unsigned int i = 0; i < 4; i++)
				   {       
				   	if(e->vn[i]->x > x_max) x_max = e->vn[i]->x;
				   	if(e->vn[i]->y > y_max) y_max = e->vn[i]->y;
				   	 if(e->vn[i]->x < x_min) x_min = e->vn[i]->x;
				   	if(e->vn[i]->y < y_min) y_min = e->vn[i]->y;
				   }
				   delta_x = (x_max - x_min)/2.; 
				   delta_y = (y_max - y_min)/2.;
				   set_delta = true;
		     }
         
         if(delta_x==0)Hermes::Exceptions::Exception("delta_x==0. ");
         if(delta_y==0)Hermes::Exceptions::Exception("delta_y==0. ");
		   
		   double u, u_1,u_2, u_3;
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
      	for(unsigned int i = 0; i < 4; i++)
         {
         	c_ref_x +=e->vn[i]->x; c_ref_y +=e->vn[i]->y; 
         	if(set_delta==false)
         	{
					if(e->vn[i]->x > x_max) x_max = e->vn[i]->x;
					if(e->vn[i]->y > y_max) y_max = e->vn[i]->y;
					 if(e->vn[i]->x < x_min) x_min = e->vn[i]->x;
					if(e->vn[i]->y < y_min) y_min = e->vn[i]->y;
         	}
         }
         c_ref_y /=4; c_ref_x/=4.; 
         if(set_delta == false)
          {
           delta_x = (x_max-x_min)/2.;
           delta_y = (y_max - y_min)/2.;
          }
          
		  for(unsigned int i = 0; i < 4; i++)
			 {
				vn = e->vn[i];
				/*u_i = find_value(e, vn->x, vn->y, 0) - u_c[e->id];
				u_dx_i = find_value(e, vn->x, vn->y, 1)- u_dx_K[e->id];
				u_dy_i = find_value(e, vn->x, vn->y, 2)- u_dy_K[e->id];	*/	
						
     
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
 		  if(alpha_i_first_order[e->id] < alpha_i_second_order[e->id])
	   		alpha_i_first_order[e->id] = alpha_i_second_order[e->id];
	   	
	   	if((alpha_i_first_order[e->id]<0.)||(alpha_i_first_order[e->id]>1.))   Hermes::Exceptions::Exception("alpha_first_order. ");
	   	if((alpha_i_second_order[e->id]<0.)||(alpha_i_second_order[e->id]>1.)) Hermes::Exceptions::Exception("alpha_sec_order. ");
      
      }
		

}


void KuzminOscillationDetector::find_centroid_values(Element* e, double &u_c)
{
    u_c = vol_average->calculate_Volume_Average(e, 0., solution, space);
  //  solution->set_active_element(e);  
   // u_c = solution->get_ref_value(e, 0., 0., 0, 0);

}

void KuzminOscillationDetector::find_centroid_derivatives(Element* e, double &u_dx_c, double &u_dy_c)
{
	//	double c_ref_x = 0.;
		//double c_ref_y = 0.; 
    solution->set_active_element(e);   
    u_dx_c = solution->get_ref_value_transformed(e, 0., 0.,  0, 1);
    u_dy_c = solution->get_ref_value_transformed(e, 0., 0.,  0, 2);
  
}

void KuzminOscillationDetector::find_second_centroid_derivatives(Element* e, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c)
{

		//double c_ref_x = 0.;
		//double c_ref_y = 0.; 
    solution->set_active_element(e);    
    u_dxx_c = solution->get_ref_value_transformed(e, 0., 0., 0, 3);
    u_dyy_c = solution->get_ref_value_transformed(e, 0., 0., 0, 4);
    u_dxy_c = solution->get_ref_value_transformed(e, 0., 0., 0, 5);    

  
}

double KuzminOscillationDetector::find_value(Element* e, double x, double y, int derivative)
{
		double c_ref_x, c_ref_y;
      solution->get_refmap()->set_active_element(e);
      solution->get_refmap()->untransform(e, x, y, c_ref_x, c_ref_y);
      return solution->get_ref_value(e, c_ref_x, c_ref_y, 0, derivative);

}

void KuzminOscillationDetector::find_vertex_values(Element* e, double vertex_values[4])
{
  double c_ref_x, c_ref_y;

    for(unsigned int j = 0; j < e->get_nvert(); j++)
    {
      solution->get_refmap()->set_active_element(e);
      solution->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_values[j] = solution->get_ref_value(e, c_ref_x, c_ref_y, 0, 0); 

      
    }
  
}

void KuzminOscillationDetector::find_vertex_derivatives(Element* e, double vertex_derivatives[4][2])
{
  double c_ref_x, c_ref_y;

    for(unsigned int j = 0; j < e->get_nvert(); j++)
    {
      solution->get_refmap()->set_active_element(e);
      solution->get_refmap()->untransform(e, e->vn[j]->x, e->vn[j]->y, c_ref_x, c_ref_y);
      vertex_derivatives[j][0] = solution->get_ref_value(e, c_ref_x, c_ref_y, 0, 1);
      vertex_derivatives[j][1] = solution->get_ref_value(e, c_ref_x, c_ref_y, 0, 2);
    }
  
}

void KuzminOscillationDetector::find_u_i_min_max_first_order(Element* e, double u_i_min[4], double u_i_max[4])
{

	Element* neigh;
  for(unsigned int j = 0; j < e->get_nvert(); j++)
  {
    NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;
    ns.set_active_edge(j);
		
    // First beginning neighbors on every edge.
    double u_c;
    ns.set_active_segment(0);
    neigh = ns.get_neighb_el();
    
    find_centroid_values(neigh, u_c);    
      if(u_i_min[j] > u_c)
        u_i_min[j] = u_c;    
      if(u_i_max[j] < u_c)
        u_i_max[j] = u_c;

    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);
    find_centroid_values(ns.get_neighb_el(), u_c);

      if(u_i_min[(j + 1) % e->get_nvert()] > u_c)
        u_i_min[(j + 1) % e->get_nvert()] = u_c;

      if(u_i_max[(j + 1) % e->get_nvert()] < u_c)
        u_i_max[(j + 1) % e->get_nvert()] = u_c;

    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert());
    ns_1.set_active_segment(0);
    
    neigh = ns_1.get_neighb_el();
        
    find_centroid_values(neigh, u_c);
      if(u_i_min[j] > u_c)
        u_i_min[j] = u_c;
      if(u_i_max[j] < u_c)
        u_i_max[j] = u_c;
  }
}

void KuzminOscillationDetector::find_alpha_i_first_order(double u_i_min[4], double u_i_max[4], double u_c, double u_i[4], double &alpha_i)
{
    alpha_i = 1;
    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
      // Sanity checks.
   /*  if(std::abs(u_i[vertex_i] - u_c) < 1E-6)
        continue;
      if(std::abs((u_i_min[vertex_i] - u_c) / u_c) > 10)
        continue;
      if(std::abs((u_i_max[vertex_i] - u_c) / u_c) > 10)
        continue;*/

		if((u_i_min[vertex_i] - u_c) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u_i_min[vertex_i] - u_c)>0. ");	
		if((u_i_max[vertex_i] - u_c) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_i_max[vertex_i] - u_c)<0"); 
		
      if(u_i[vertex_i] < u_c)
      {
        if((u_i_min[vertex_i] - u_c) / (u_i[vertex_i] - u_c) < alpha_i)
          alpha_i = (u_i_min[vertex_i] - u_c) / (u_i[vertex_i] - u_c);
      }
      else if(u_i[vertex_i] > u_c)
      {
        if((u_i_max[vertex_i] - u_c) / (u_i[vertex_i] - u_c) < alpha_i)
          alpha_i = (u_i_max[vertex_i] - u_c) / (u_i[vertex_i] - u_c);
      }
    }
  
}


void KuzminOscillationDetector::find_u_i_min_max_second_order(Element* e, double u_d_i_min[4][2], double u_d_i_max[4][2])
{
  for(unsigned int j = 0; j < e->get_nvert(); j++)
  {
    NeighborSearch<double> ns(e, mesh);
    if(e->en[j]->bnd)
      continue;

    ns.set_active_edge(j);

    // First beginning neighbors on every edge.
    double u_dx_c, u_dy_c;
    ns.set_active_segment(0);
    find_centroid_derivatives(ns.get_neighb_el(), u_dx_c, u_dy_c);
    
    ////-----------
    
/*    double delta_x = get_delta_x(ns.get_neighb_el());
    double delta_y = get_delta_y(ns.get_neighb_el());
    
    u_dx_c*=delta_x; u_dy_c*=delta_y;*/
    ///------------

      if(u_d_i_min[j][0] > u_dx_c)
        u_d_i_min[j][0] = u_dx_c;
      if(u_d_i_min[j][1] > u_dy_c)
        u_d_i_min[j][1] = u_dy_c;
   

      if(u_d_i_max[j][0] < u_dx_c)
        u_d_i_max[j][0] = u_dx_c;
      if(u_d_i_max[j][1] < u_dy_c)
        u_d_i_max[j][1] = u_dy_c;
   
    // Second end neighbors on every edge.
    ns.set_active_segment(ns.get_num_neighbors() - 1);
    find_centroid_derivatives(ns.get_neighb_el(), u_dx_c, u_dy_c);
    
        ////-----------
    
 /*   delta_x = get_delta_x(ns.get_neighb_el());
    delta_y = get_delta_y(ns.get_neighb_el());
    
    u_dx_c*=delta_x; u_dy_c*=delta_y;*/
    ///------------

      if(u_d_i_min[(j + 1) % e->get_nvert()][0] > u_dx_c)
        u_d_i_min[(j + 1) % e->get_nvert()][0] = u_dx_c;
      if(u_d_i_min[(j + 1) % e->get_nvert()][1] > u_dy_c)
        u_d_i_min[(j + 1) % e->get_nvert()][1] = u_dy_c;

      if(u_d_i_max[(j + 1) % e->get_nvert()][0] < u_dx_c)
        u_d_i_max[(j + 1) % e->get_nvert()][0] = u_dx_c;
      if(u_d_i_max[(j + 1) % e->get_nvert()][1] < u_dy_c)
        u_d_i_max[(j + 1) % e->get_nvert()][1] = u_dy_c;


    // Now the hard part, neighbors' neighbors.
    /// \todo This is where it fails for triangles, where it is much more complicated to look for elements sharing a vertex.
    ns.set_active_segment(0);
    NeighborSearch<double> ns_1(ns.get_neighb_el(), mesh);
    if(ns.get_neighb_el()->en[(ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert()]->bnd)
      continue;
    ns_1.set_active_edge((ns.get_neighbor_edge().local_num_of_edge + 1) % ns.get_neighb_el()->get_nvert());
    ns_1.set_active_segment(0);
    find_centroid_derivatives(ns_1.get_neighb_el(), u_dx_c, u_dy_c);
    
            ////-----------
    
  /*  delta_x = get_delta_x(ns_1.get_neighb_el());
    delta_y = get_delta_y(ns_1.get_neighb_el());
    
    u_dx_c*=delta_x; u_dy_c*=delta_y;*/
    ///------------

      if(u_d_i_min[j][0] > u_dx_c)
        u_d_i_min[j][0] = u_dx_c;
      if(u_d_i_min[j][1] > u_dy_c)
        u_d_i_min[j][1] = u_dy_c;

      if(u_d_i_max[j][0] < u_dx_c)
        u_d_i_max[j][0] = u_dx_c;
      if(u_d_i_max[j][1] < u_dy_c)
        u_d_i_max[j][1] = u_dy_c;
    
  }
}

void KuzminOscillationDetector::find_alpha_i_second_order(double u_d_i_min[4][2], double u_d_i_max[4][2], double &u_dx_c, double &u_dy_c, double u_d_i[4][2], double &alpha_i)
{
    alpha_i = 1;
    for(unsigned int vertex_i = 0; vertex_i < 4; vertex_i++)
    {  

    
      // Sanity checks.
   /*  if(std::abs(u_dx_c - u_d_i[vertex_i][0]) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[vertex_i][0] - u_dx_c) / u_dx_c) > 10)
        continue;
      if(std::abs((u_d_i_max[vertex_i][0] - u_dx_c) / u_dx_c) > 10)
        continue;

      if(std::abs(u_dy_c - u_d_i[vertex_i][1]) < 1E-5)
        continue;
      if(std::abs((u_d_i_min[vertex_i][1] - u_dy_c) / u_dy_c) > 10)
        continue;
      if(std::abs((u_d_i_max[vertex_i][1] - u_dy_c) / u_dy_c) > 10)
        continue;*/

		if((u_d_i_min[vertex_i][0] - u_dx_c) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u__d_i_min[vertex_i][0] - u__dx_c)>0. ");	
		if((u_d_i_max[vertex_i][0] - u_dx_c) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_di_max[vertex_i][0] - u_dx_c)<0"); 
		if((u_d_i_min[vertex_i][1] - u_dy_c) >0) Hermes::Exceptions::Exception(" KuzminOscillationDetector::(u__d_i_min[vertex_i][1] - u__dy_c)>0. ");	
		if((u_d_i_max[vertex_i][1] - u_dy_c) <0) Hermes::Mixins::Loggable::Static::info("KuzminOscillationDetector:: (u_di_max[vertex_i][1] - u_dy_c)<0"); 

      // dx.
      if(u_d_i[vertex_i][0] < u_dx_c)
      {
        if((u_d_i_min[vertex_i][0] - u_dx_c) / (u_d_i[vertex_i][0] - u_dx_c) < alpha_i)
          alpha_i = (u_d_i_min[vertex_i][0] - u_dx_c) / (u_d_i[vertex_i][0] - u_dx_c);
      }
      else if(u_d_i[vertex_i][0] > u_dx_c)
      {
        if((u_d_i_max[vertex_i][0] - u_dx_c) / (u_d_i[vertex_i][0] - u_dx_c) < alpha_i)
          alpha_i = (u_d_i_max[vertex_i][0] - u_dx_c) / (u_d_i[vertex_i][0] - u_dx_c);
      }
      // dy.
      if(u_d_i[vertex_i][1] < u_dy_c)
      {
        if((u_d_i_min[vertex_i][1] - u_dy_c) / (u_d_i[vertex_i][1] - u_dy_c) < alpha_i)
          alpha_i = (u_d_i_min[vertex_i][1] - u_dy_c) / (u_d_i[vertex_i][1] - u_dy_c);
      }
      else if(u_d_i[vertex_i][1] > u_dy_c)
      {
        if((u_d_i_max[vertex_i][1] - u_dy_c) / (u_d_i[vertex_i][1] - u_dy_c) < alpha_i)
          alpha_i = (u_d_i_max[vertex_i][1] - u_dy_c) / (u_d_i[vertex_i] [1]- u_dy_c);
      }
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
  
  
  
  
  
  SlopeLimiterSolution::SlopeLimiterSolution(Solution<double>* ref_sln, Space<double>* space):Solution<double>(space->get_mesh()),ref_sln(ref_sln), space(space) 
   { 
  detector = new KuzminOscillationDetector(space, ref_sln); 
  };
  
 
    
    
 void SlopeLimiterSolution::limit_solution_according_to_detector()
    {
     int o;			
		this->increasePointerDataCounter();
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

        // Hcurl and Hdiv: actual order of functions is one higher than element order
        if((space->get_shapeset())->get_num_components() == 2) o++;

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
      Mesh* ref_mesh = new Mesh;
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
      double u_c, u_dx_c, u_dy_c, u_dxx_c, u_dxy_c, u_dyy_c;   
           
      detector->compute_all();      
     refmap->set_quad_2d(quad);
        
      for_all_active_elements(e, this->mesh)
      {
     		 shapeset_taylor->init_data(e, space);
          PrecalcShapeset *pss_taylor = new PrecalcShapeset(shapeset_taylor);
          pss_taylor->set_quad_2d(quad);
          
	 		 ref_space->set_shapeset(shapeset_taylor);
	     	 ref_space->assign_dofs();

      
		  // detector->get_alpha_all_orders(e,alpha, u_c, u_dx_c, u_dy_c, u_dxx_c, u_dxy_c, u_dyy_c);
			//lokale Taylorbasisfunktionen ohne delta_x, delta_y!!     	 
			
	 		alpha[0] = detector->get_alpha_first_order()[e->id];
			 alpha[1] = detector->get_alpha_second_order()[e->id];  
   		/*		u_c = detector->get_u_c(e->id);
 			u_dx_c= detector->get_u_dx_K(e->id)/shapeset->get_delta_x(); 
  			u_dy_c= detector->get_u_dy_K(e->id)/shapeset->get_delta_y(); 
  			u_dxx_c= detector->get_u_dxx_K(e->id)/(shapeset->get_delta_x()*shapeset->get_delta_x());
  			u_dyy_c= detector->get_u_dyy_K(e->id)/(shapeset->get_delta_y()*shapeset->get_delta_y());  
  			u_dxy_c= detector->get_u_dxy_K(e->id)/(shapeset->get_delta_x()*shapeset->get_delta_y());*/
  			
  			detector->find_centroid_values(e, u_c);
  			detector->find_centroid_derivatives(e, u_dx_c, u_dy_c);
 			detector->find_second_centroid_derivatives(e, u_dxx_c, u_dxy_c, u_dyy_c);

			
			//alpha[0] = 0.; 
			//alpha[1] = 0.;
			
        this->mode = e->get_mode();
        o = elem_orders[e->id];
        int np = quad->get_num_points(o, e->get_mode());

        AsmList<double> al;
        ref_space->get_element_assembly_list(e, &al);
        pss_taylor->set_active_element(e);             						
						
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
			  
				pss_taylor->set_active_shape(al.get_idx()[k]);
				pss_taylor->set_quad_order(o, H2D_FN_VAL);      

				double coef = al.get_coef()[k] * u_new;  
				double3* pts  = quad->get_points(o, e->get_mode());				



				//double* shape = pss_taylor->get_fn_values(l);
				for (int i = 0; i < np; i++)
				{
					//val[i] += shape[i] * coef;				
				val[i] += pts[i][2]*shapeset_taylor->get_value(0, al.get_idx()[k], x[i], y[i], 0, e->get_mode())* coef;					
					//val[i] += pts[i][2]*shapeset_taylor->get_value(0, al.get_idx()[k], pts[i][0], pts[i][1], 0, e->get_mode())* coef;	   
				}        
          }
          mono += np;

          // solve for the monomial coefficients
          if(mono_lu.mat[this->mode][o] == NULL)
            mono_lu.mat[this->mode][o] = calc_mono_matrix(o, mono_lu.perm[this->mode][o]);
          lubksb(mono_lu.mat[this->mode][o], np, mono_lu.perm[this->mode][o], val);
        }
        delete pss_taylor;
      }

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL.\n");
      init_dxdy_buffer();
      this->element = NULL;
       
      
        delete [] alpha;
        delete ref_mesh;
        delete ref_space;
        delete shapeset_taylor;
       
    }

