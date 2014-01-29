#include "interface.h"

/////////////------------WEAKFORM-------------


	EulerInterface::EulerInterface(double kappa, MeshSharedPtr mesh, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,NumericalFlux* num_flux, EulerFluxes* euler_fluxes,
RiemannInvariants* riemann_invariants, int num_of_equations): WeakForm<double>(num_of_equations), num_flux(num_flux),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mesh(mesh)
	{
	
		for(int form_i = 0; form_i < 4; form_i++)
		{				
				  add_vector_form_DG(new EulerInterface::EulerEquationsVectorFormFlux(form_i, kappa,num_flux));
			for(int form_j = 0; form_j<4; form_j++)
				add_matrix_form_DG(new EulerInterface::EulerEquationsBilinearFormFlux(form_i, form_j, kappa, num_flux));
		}

   this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

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
      double w_L[4], w_R[4],jump_u[4], mid_u[4];
		double nx, ny, tx, ty;

      for (int point_i = 0; point_i < n; point_i++)
      {

		nx = e->nx[point_i];
		ny = e->ny[point_i];
		tx = e->tx[point_i];
		ty = e->ty[point_i];
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];	

		
        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, nx, ny) * v->val[point_i];


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
        double result = 0.;
      double w_L[4], w_R[4],w_mid[4], w_mean[4];
	double * A_n = new double[4];
double nx, ny, tx, ty;

      for (int point_i = 0; point_i < n; point_i++)
      {
		double jump_v = (v->fn_central == NULL ? -v->val_neighbor[point_i] : v->val[point_i]);
		double jump_u =(u->fn_central == NULL ? -u->val_neighbor[point_i] :u->val[point_i]);
		double mid_u = 0.5* (u->fn_central == NULL ? u->val_neighbor[point_i] :u->val[point_i]);

		nx = e->nx[point_i];
		ny = e->ny[point_i];
		tx = e->tx[point_i];
		ty = e->ty[point_i];

        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];

		w_mid[0] = 0.5*(w_L[0]+w_R[0]);
		w_mid[1] = 0.5*(w_L[1]+w_R[1]);
		w_mid[2] = 0.5*(w_L[2]+w_R[2]);
		w_mid[3] = 0.5*(w_L[3]+w_R[3]);


	if(dynamic_cast<LaxFriedrichsNumericalFlux* >(this->num_flux) != NULL)
	{

	double s_right = std::fabs(((nx * w_R[1]) + (ny * w_R[2])) / w_R[0]) + QuantityCalculator::calc_sound_speed(w_R[0], w_R[1], w_R[2], w_R[3], this->kappa);
	double s_left = std::fabs(((nx * w_L[1]) + (ny * w_L[2])) / w_L[0]) + QuantityCalculator::calc_sound_speed(w_L[0], w_L[1], w_L[2], w_L[3], this->kappa);

	if(this->i==this->j) result += wt[point_i]*jump_v*jump_u*std::max(s_left, s_right)*0.5;


		if(u->fn_central != NULL){
			result += wt[point_i]*jump_v * mid_u * 
				( (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_L[0], w_L[1], w_L[2], w_L[3],0,this->i,this->j) 
				  * nx+        
				 (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_L[0], w_L[1], w_L[2], w_L[3],1,this->i,this->j) 
				  * ny);

		}else{
			result += wt[point_i]*jump_v * mid_u * 
				( (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_R[0], w_R[1], w_R[2], w_R[3],0,this->i,this->j) 
				  * nx+        
				 (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_R[0], w_R[1], w_R[2], w_R[3],1,this->i,this->j) 
				  * ny	);
		}

	}else if(dynamic_cast<VijayasundaramNumericalFlux* >(this->num_flux) != NULL)
	{		 
		if(u->fn_central != NULL)	
			Boundary_helpers::calculate_P_plus(w_L[0], w_L[1], w_L[2], w_L[3], nx,ny , w_R[0], w_R[1], w_R[2], w_R[3], kappa, this->i,A_n);
		else
			Boundary_helpers::calculate_P_minus(w_L[0], w_L[1], w_L[2], w_L[3], nx,ny , w_R[0], w_R[1], w_R[2], w_R[3], kappa, this->i,A_n);

		result += wt[point_i]*jump_v*mid_u*2.* A_n[this->j];
			  
	}else if(dynamic_cast<ApproxRoeNumericalFlux* >(this->num_flux) != NULL) {

		Boundary_helpers::calculate_A_n(w_L[0], w_L[1], w_L[2], w_L[3], nx,ny , w_R[0], w_R[1], w_R[2], w_R[3], kappa, this->i,A_n);
		result += wt[point_i]*jump_v*jump_u*0.5* A_n[this->j];

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
  return new EulerInterface::EulerEquationsBilinearFormFlux(this->i,this->j, this->kappa, this->num_flux);
}



