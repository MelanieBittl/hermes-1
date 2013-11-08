#include "interface.h"

/////////////------------WEAKFORM-------------


	EulerInterface::EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy,NumericalFlux* num_flux, int num_of_equations): WeakForm<double>(num_of_equations), num_flux(num_flux)
	{

		for(int form_i = 0; form_i < 4; form_i++)
		{
				EulerEquationsVectorFormFlux* formDG = new EulerEquationsVectorFormFlux(form_i, kappa,num_flux);
				  add_vector_form_DG(formDG);
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
double EulerEquationsVectorFormFlux::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, 
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
