#include "interface.h"

/////////////------------WEAKFORM-------------


	EulerInterface::EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,NumericalFlux* num_flux, int num_of_equations): WeakForm<double>(num_of_equations), num_flux(num_flux),euler_fluxes(new EulerFluxes(kappa))
	{
	
		for(int form_i = 0; form_i < 4; form_i++)
		{				
				  add_vector_form_DG(new EulerInterface::EulerEquationsVectorFormFlux(form_i, kappa,num_flux));
			for(int form_j = 0; form_j<4; form_j++)
				add_matrix_form_DG(new EulerInterface::EulerEquationsBilinearFormFlux(form_i, form_j, kappa));
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
      double w_L[4], w_R[4];
      for (int point_i = 0; point_i < n; point_i++)
      {
        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];

        result += wt[point_i] * this->num_flux->numerical_flux_i(this->i, w_L, w_R, e->nx[point_i], e->ny[point_i]) * v->val[point_i];
      }

      return -result;
    }


    Ord EulerInterface::EulerEquationsVectorFormFlux::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const
{  
Ord result = Ord(10); 
  return result;
}





//----------DG-Matrix-Form------------
/*
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
	//Real v_x =Real(1.);
 	//Real v_y =Real(1.);
		Real v_x = (0.5- e->y[i]); 
 		Real v_y = (e->x[i]-0.5) ; 
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(Real(0), u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Real(0), a_dot_n) * jump_v;
      
  }
  return result*theta;
}*/

double EulerInterface::EulerEquationsBilinearFormFlux::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
        double result = 0.;
      double w_L[4], w_R[4];
      for (int point_i = 0; point_i < n; point_i++)
      {
		double jump_v = (v->fn_central == NULL ? -v->val_neighbor[point_i] : v->val[point_i]);
		double jump_u =(u->fn_central == NULL ? -u->val_neighbor[point_i] :u->val[point_i]);

        w_L[0] = ext[0]->val[point_i];
        w_L[1] = ext[1]->val[point_i];
        w_L[2] = ext[2]->val[point_i];
        w_L[3] = ext[3]->val[point_i];

        w_R[0] = ext[0]->val_neighbor[point_i];
        w_R[1] = ext[1]->val_neighbor[point_i];
        w_R[2] = ext[2]->val_neighbor[point_i];
        w_R[3] = ext[3]->val_neighbor[point_i];


double s_right = std::abs(((e->nx[point_i] * w_R[1]) + (e->ny[point_i] * w_R[2])) / w_R[0]) + QuantityCalculator::calc_sound_speed(w_R[0], w_R[1], w_R[2], w_R[3], this->kappa);
double s_left = std::abs(((e->nx[point_i] * w_L[1]) + (e->ny[point_i] * w_L[2])) / w_L[0]) + QuantityCalculator::calc_sound_speed(w_L[0], w_L[1], w_L[2], w_L[3], this->kappa);

    if(u->fn_central != NULL)
        result += wt[point_i]*0.5*jump_v *( u->val[point_i] * 
        ( (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_L[0], w_L[1], w_L[2], w_L[3],0,this->i,this->j) 
          * e->nx[point_i]+        
         (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_L[0], w_L[1], w_L[2], w_L[3],1,this->i,this->j) 
          * e->ny[point_i]));
	else
 		result += wt[point_i]*0.5*jump_v *(u->val_neighbor[point_i]*
			  ( (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_R[0], w_R[1], w_R[2], w_R[3],0,this->i,this->j) 
					  * e->nx[point_i]+        
					 (static_cast<EulerInterface*>(wf))->euler_fluxes->A(w_R[0], w_R[1], w_R[2], w_R[3],1,this->i,this->j) 
					  * e->ny[point_i]) );
		
		result += wt[point_i]*jump_v*jump_u*std::max(s_left, s_right); 
		

      
      }
return -result;
}

Ord EulerInterface::EulerEquationsBilinearFormFlux::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
 Ord result = Ord(10); 
  return result;
}

MatrixFormDG<double>* EulerInterface::EulerEquationsBilinearFormFlux::clone() const
{
  return new EulerInterface::EulerEquationsBilinearFormFlux(this->i,this->j, this->kappa);
}



