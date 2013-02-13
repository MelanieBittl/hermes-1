#include "util.h" 
//derivative 0 = u, derivative 1 = du/dx, derivative 2 = du/dy

double VolumeAverage::calculate_Volume_Average(Element* e,int derivative, Solution<double>* sln_1,
Space<double>* ref_space, Solution<double>* sln_2)
{
	if(e==NULL) Hermes::Exceptions::Exception("Volume_Average: Element* e ==NULL");
	if(sln_1==NULL) Hermes::Exceptions::Exception("Volume_Average: sln_1 ==NULL");	
	if(derivative<0) Hermes::Exceptions::Exception("Volume_Average: derivative<0?");	
	if(ref_space==NULL) Hermes::Exceptions::Exception("Volume_Average: space=NULL?");	
	
	double result = 0;
	int ndof = ref_space->get_num_dofs();
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
			if(derivative==0)
			{
				if(sln_2==NULL)
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->val[0];	
				else 
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->val[0]*(sln_2->get_pt_value(x[i],y[i]))->val[0];
			}else if(derivative==1)
			{	if(sln_2==NULL)
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->dx[0];	
				else 
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->dx[0]*(sln_2->get_pt_value(x[i],y[i]))->val[0];				
			}else if(derivative==2)
			{
				if(sln_2==NULL)
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->dy[0];	
				else 
					result += pts[i][2]*(sln_1->get_pt_value(x[i],y[i]))->dy[0]*(sln_2->get_pt_value(x[i],y[i]))->val[0];	
			}else{
				Hermes::Exceptions::Exception("Volume_Average: for HERMES_EXACT no higher derivatives ");			
			}
		}
		delete refmap;
	
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

	result /= std::abs(e->get_area()) ;
	result *=jac;
	return result;
}


