#include "util.h" 
//derivative 0 = u, derivative 1 = du/dx, derivative 2 = du/dy

double VolumeAverage::calculate_Volume_Average(Element* e,int derivative, MeshFunctionSharedPtr<double>  sln_1_fct,
SpaceSharedPtr<double>  ref_space, MeshFunctionSharedPtr<double>  sln_2_fct)
{
	if(e==NULL) Hermes::Exceptions::Exception("Volume_Average: Element* e ==NULL");
	if(sln_1_fct==NULL) Hermes::Exceptions::Exception("Volume_Average: sln_1 ==NULL");	
	if(derivative<0) Hermes::Exceptions::Exception("Volume_Average: derivative<0?");	
	if(ref_space==NULL) Hermes::Exceptions::Exception("Volume_Average: space=NULL?");	

Solution<double>* sln_1 = (dynamic_cast<Solution<double>*>(sln_1_fct.get()));
Solution<double>* sln_2 =NULL;
if(sln_2_fct!=NULL)
 sln_2 = (dynamic_cast<Solution<double>*>(sln_2_fct.get()));

	
	double result = 0;
	//int ndof = ref_space->get_num_dofs();
	int o = ref_space->get_element_order(e->id);
  o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o)); 
  for (unsigned int k = 0; k < e->get_nvert(); k++)
  {
	 int eo = ref_space->get_edge_order(e, k);
	 if(eo > o) o = eo;
  }

		RefMap* refmap = new RefMap();
		refmap->set_quad_2d(&g_quad_2d_std);
		refmap->set_active_element(e);		
      double jac = refmap->get_const_jacobian();

	if(sln_1->get_type() ==HERMES_EXACT)
	{		
		double* x =  refmap->get_phys_x(o);
		double* y =  refmap->get_phys_y(o);
		Quad2D* quad = refmap->get_quad_2d();
		int np = quad->get_num_points(o, e->get_mode());
		double3* pts  = quad->get_points(o, e->get_mode());
		for (int i = 0; i < np; i++)
		{
			Func<double>* s1 = sln_1->get_pt_value(x[i],y[i]);
			Func<double>* s2 = NULL;
			if(sln_2!=NULL) s2 = sln_2->get_pt_value(x[i],y[i]);
			if(derivative==0)
			{
				if(s2==NULL)
					result += pts[i][2]*s1->val[0];	
				else 
					result += pts[i][2]*s1->val[0]*s2->val[0];
			}else if(derivative==1)
			{	if(s2==NULL)
					result += pts[i][2]*s1->dx[0];	
				else 
					result += pts[i][2]*s1->dx[0]*s2->val[0];				
			}else if(derivative==2)
			{
				if(s2==NULL)
					result += pts[i][2]*s1->dy[0];	
				else 
					result += pts[i][2]*s1->dy[0]*s2->val[0];	
			}else{
				Hermes::Exceptions::Exception("Volume_Average: for HERMES_EXACT no higher derivatives ");			
			}
			s1->free_fn(); delete s1;
			if(s2!=NULL)
			{
				 s2->free_fn();
				delete s2;	 
			}
		}
		
	
	}else{	
		Quad2D* quad = &g_quad_2d_std;
		int np = quad->get_num_points(o, e->get_mode());
		double3* pts  = quad->get_points(o, e->get_mode());

		for (int i = 0; i < np; i++)
		{
			if(sln_2==NULL)
				result += pts[i][2]*sln_1->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative);	
			else 
				result += pts[i][2]*sln_1->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative)*sln_2->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative);
		}
	}	
delete refmap;
	result /= std::abs(e->get_area()) ;
	result *=jac;

	return result;
}

 double VolumeAverage::calculate_Volume_Average(Element* e,int derivative,  Solution<double>* sln_1, SpaceSharedPtr<double> ref_space,  Solution<double>* sln_2)
{
	if(e==NULL) Hermes::Exceptions::Exception("Volume_Average: Element* e ==NULL");
	if(sln_1==NULL) Hermes::Exceptions::Exception("Volume_Average: sln_1 ==NULL");	
	if(derivative<0) Hermes::Exceptions::Exception("Volume_Average: derivative<0?");	
	if(ref_space==NULL) Hermes::Exceptions::Exception("Volume_Average: space=NULL?");	

	
	double result = 0;
	//int ndof = ref_space->get_num_dofs();
	int o = ref_space->get_element_order(e->id);
  o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o)); 
  for (unsigned int k = 0; k < e->get_nvert(); k++)
  {
	 int eo = ref_space->get_edge_order(e, k);
	 if(eo > o) o = eo;
  }

		RefMap* refmap = new RefMap();
		refmap->set_quad_2d(&g_quad_2d_std);
		refmap->set_active_element(e);		
      double jac = refmap->get_const_jacobian();

	if(sln_1->get_type() ==HERMES_EXACT)
	{		
		double* x =  refmap->get_phys_x(o);
		double* y =  refmap->get_phys_y(o);
		Quad2D* quad = refmap->get_quad_2d();
		int np = quad->get_num_points(o, e->get_mode());
		double3* pts  = quad->get_points(o, e->get_mode());
		for (int i = 0; i < np; i++)
		{
			Func<double>* s1 = sln_1->get_pt_value(x[i],y[i]);
			Func<double>* s2 = NULL;
			if(sln_2!=NULL) s2 = sln_2->get_pt_value(x[i],y[i]);
			if(derivative==0)
			{
				if(s2==NULL)
					result += pts[i][2]*s1->val[0];	
				else 
					result += pts[i][2]*s1->val[0]*s2->val[0];
			}else if(derivative==1)
			{	if(s2==NULL)
					result += pts[i][2]*s1->dx[0];	
				else 
					result += pts[i][2]*s1->dx[0]*s2->val[0];				
			}else if(derivative==2)
			{
				if(s2==NULL)
					result += pts[i][2]*s1->dy[0];	
				else 
					result += pts[i][2]*s1->dy[0]*s2->val[0];	
			}else{
				Hermes::Exceptions::Exception("Volume_Average: for HERMES_EXACT no higher derivatives ");			
			}
			s1->free_fn(); delete s1;
			if(s2!=NULL)
			{
				 s2->free_fn();
				delete s2;	 
			}
		}
		
	
	}else{	
		Quad2D* quad = &g_quad_2d_std;
		int np = quad->get_num_points(o, e->get_mode());
		double3* pts  = quad->get_points(o, e->get_mode());

		for (int i = 0; i < np; i++)
		{
			if(sln_2==NULL)
				result += pts[i][2]*sln_1->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative);	
			else 
				result += pts[i][2]*sln_1->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative)*sln_2->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, derivative);
		}
	}	
delete refmap;
	result /= std::abs(e->get_area()) ;
	result *=jac;

	return result;
}


 double VolumeAverage::calculate_norm(Element* e,SpaceSharedPtr<double>  space, MeshFunctionSharedPtr<double>  sln_fct)
 {
 
 	if(e==NULL) Hermes::Exceptions::Exception("Volume_Average:calculate_norm Element* e ==NULL");
	if(sln_fct==NULL) Hermes::Exceptions::Exception("Volume_Average:calculate_norm sln_1 ==NULL");
	if(space==NULL) Hermes::Exceptions::Exception("Volume_Average:calculate_norm space=NULL?");		

	Solution<double>* sln = (dynamic_cast<Solution<double>*>(sln_fct.get()));
	double result = 0;
	int o = space->get_element_order(e->id);
  o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o)); 

		RefMap* refmap = new RefMap;
		refmap->set_quad_2d(&g_quad_2d_std);
		refmap->set_active_element(e);		
      double jac = refmap->get_const_jacobian();
	
		Quad2D* quad = &g_quad_2d_std;
		int np = quad->get_num_points(o, e->get_mode());
		double3* pts  = quad->get_points(o, e->get_mode());

		for (int i = 0; i < np; i++)
				result += pts[i][2]*sln->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, 0)*sln->get_ref_value_transformed(e, pts[i][0], pts[i][1], 0, 0);
		
		
	result *=jac;
	delete refmap;
	return result;
 
 
 }
