#include "interface.h"

/////////////------------WEAKFORM-------------

double calculate_a_dot_v(double x, double y, double vx, double vy)
{
 
 return  x*vx + y*vy;
}


double upwind_flux(double u_cent, double u_neib, double a_dot_n) 
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}




	EulerInterface::EulerInterface(double gamma, MeshSharedPtr mesh, MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p,NumericalFlux* num_flux, EulerFluxes* euler_fluxes,
RiemannInvariants* riemann_invariants, int num_of_equations): WeakForm<double>(num_of_equations), num_flux(num_flux),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mesh(mesh)
	{
	//gas/particle phase
		for(int form_i = 0; form_i < 5; form_i++)
		{				
			add_vector_form_DG(new EulerInterface::EulerEquationsVectorFormFlux(form_i, gamma,num_flux));
			for(int form_j = 0; form_j<5; form_j++)
				add_matrix_form_DG(new EulerInterface::EulerEquationsBilinearFormFlux(form_i, form_j, gamma, num_flux));
		}




    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density_g, prev_density_vel_x_g, prev_density_vel_y_g, prev_energy_g,
				prev_density_p));

	};


	EulerInterface ::~EulerInterface ()
	{
	};
	
	    WeakForm<double>* EulerInterface::clone() const
    {
      const_cast<EulerInterface*>(this)->warned_nonOverride = false;
      return new EulerInterface(*this);
    }




/////////////////----------------------Euler Vectorform--------------------------------
double EulerInterface::EulerEquationsVectorFormFlux::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
      Func<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double result = 0.;
      double w_L[4], w_R[4];
		double nx, ny, tx, ty;

      for (int point_i = 0; point_i < n; point_i++)
      {

		nx = e->nx[point_i];
		ny = e->ny[point_i];
		tx = e->tx[point_i];
		ty = e->ty[point_i];	

		if(entry_i<4) //gas phase
		{
				w_L[0] = ext[0]->val[point_i];
				w_L[1] = ext[1]->val[point_i];
				w_L[2] = ext[2]->val[point_i];
				w_L[3] = ext[3]->val[point_i];

				w_R[0] = ext[0]->val_neighbor[point_i];
				w_R[1] = ext[1]->val_neighbor[point_i];
				w_R[2] = ext[2]->val_neighbor[point_i];
				w_R[3] = ext[3]->val_neighbor[point_i];
				result += wt[point_i] * this->num_flux->numerical_flux_i(entry_i, w_L, w_R, nx, ny) * v->val[point_i];
		}else //particle
		{

			w_L[0] = ext[4]->val[point_i];
			w_R[0] = ext[4]->val_neighbor[point_i];

		double v_x = 0.5*(ext[1]->val[point_i]/ext[0]->val[point_i]+ext[1]->val_neighbor[point_i]/ext[0]->val_neighbor[point_i]); 
 		double v_y = 0.5*(ext[2]->val[point_i]/ext[0]->val[point_i]+ext[2]->val_neighbor[point_i]/ext[0]->val_neighbor[point_i]); 

   double a_dot_n = calculate_a_dot_v(v_x, v_y, e->nx[point_i], e->ny[point_i]);

    double jump_v =  v->val[point_i];

      result += wt[point_i] * upwind_flux(w_L[0], w_R[0], a_dot_n) * jump_v;


		}


     }

      return (-result);
    }


    Ord EulerInterface::EulerEquationsVectorFormFlux::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const
{  
Ord result = Ord(10); 
  return result;
}





//----------DG-Matrix-Form------------


double EulerInterface::EulerEquationsBilinearFormFlux::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{

if((entry_i<4)&&(entry_j>=4)) return 0;
if((entry_i>=4)&&(entry_j<4)) return 0;
        double result = 0.;
      double w_L[4], w_R[4],w_mid[4], w_mean[4];
	double * A_n = new double[4];
double nx, ny, tx, ty, s_left, s_right; bool particle = false;

  for (int point_i = 0; point_i < n; point_i++)
      {
		double jump_v = (v->fn_central == NULL ? -v->val_neighbor[point_i] : v->val[point_i]);
		double jump_u =(u->fn_central == NULL ? -u->val_neighbor[point_i] :u->val[point_i]);
		double mid_u = 0.5* (u->fn_central == NULL ? u->val_neighbor[point_i] :u->val[point_i]);

		nx = e->nx[point_i];
		ny = e->ny[point_i];
		tx = e->tx[point_i];
		ty = e->ty[point_i];
		if((entry_i<4)&&(entry_j<4)) //gas phase
  		{   particle = false;  
			w_L[0] = ext[0]->val[point_i];
		    w_L[1] = ext[1]->val[point_i];
		    w_L[2] = ext[2]->val[point_i];
		    w_L[3] = ext[3]->val[point_i];

		    w_R[0] = ext[0]->val_neighbor[point_i];
		    w_R[1] = ext[1]->val_neighbor[point_i];
		    w_R[2] = ext[2]->val_neighbor[point_i];
		    w_R[3] = ext[3]->val_neighbor[point_i];
		}else{
			particle = true;
	 		w_L[0] = ext[4]->val[point_i];
		    w_R[0] = ext[4]->val_neighbor[point_i];

		}

		w_mid[0] = 0.5*(w_L[0]+w_R[0]);
		w_mid[1] = 0.5*(w_L[1]+w_R[1]);
		w_mid[2] = 0.5*(w_L[2]+w_R[2]);
		w_mid[3] = 0.5*(w_L[3]+w_R[3]);


	if(dynamic_cast<LaxFriedrichsNumericalFlux* >(this->num_flux) != NULL)
	{
		if(particle)
		{
				double v_x = 0.5*(ext[1]->val[point_i]/ext[0]->val[point_i]+ext[1]->val_neighbor[point_i]/ext[0]->val_neighbor[point_i]); 
 				double v_y = 0.5*(ext[2]->val[point_i]/ext[0]->val[point_i]+ext[2]->val_neighbor[point_i]/ext[0]->val_neighbor[point_i]); 

			   double a_dot_n = calculate_a_dot_v(v_x, v_y, e->nx[point_i], e->ny[point_i]);
			 
				if(u->fn_central == NULL)
				  result += wt[point_i] * upwind_flux(0, u->val_neighbor[point_i], a_dot_n) * jump_v;
				else
				  result += wt[point_i] * upwind_flux(u->val[point_i],0, a_dot_n) * jump_v;
		}else{
			s_right = std::fabs(((nx * w_R[1]) + (ny * w_R[2])) / w_R[0]) + QuantityCalculator::calc_sound_speed(w_R[0], w_R[1], w_R[2], w_R[3], this->gamma);
			s_left = std::fabs(((nx * w_L[1]) + (ny * w_L[2])) / w_L[0]) + QuantityCalculator::calc_sound_speed(w_L[0], w_L[1], w_L[2], w_L[3], this->gamma);
		

			if(entry_i==entry_j) result += wt[point_i]*jump_v*jump_u*std::max(s_left, s_right)*0.5;


			if(u->fn_central != NULL){
				result += wt[point_i]*jump_v * mid_u * 
					( (static_cast<EulerInterface*>(wf))->euler_fluxes->A_g(w_L[0], w_L[1], w_L[2], w_L[3],0,entry_i,entry_j) 
					  * nx+        
					 (static_cast<EulerInterface*>(wf))->euler_fluxes->A_g(w_L[0], w_L[1], w_L[2], w_L[3],1,entry_i,entry_j) 
					  * ny);
			}else{

				result += wt[point_i]*jump_v * mid_u * 
					( (static_cast<EulerInterface*>(wf))->euler_fluxes->A_g(w_R[0], w_R[1], w_R[2], w_R[3],0,entry_i,entry_j) 
					  * nx+        
					 (static_cast<EulerInterface*>(wf))->euler_fluxes->A_g(w_R[0], w_R[1], w_R[2], w_R[3],1,entry_i,entry_j) 
					  * ny	);
			}
		}

	

	}

}

delete [] A_n;
return -result;
}

Ord EulerInterface::EulerEquationsBilinearFormFlux::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
 Ord result = Ord(10); 
  return result;
}

MatrixFormDG<double>* EulerInterface::EulerEquationsBilinearFormFlux::clone() const
{
  return new EulerInterface::EulerEquationsBilinearFormFlux(entry_i,entry_j, this->gamma, this->num_flux);
}



