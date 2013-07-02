#include "interface.h"


    double EulerInterfaceBilinearForm::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, 
      DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double>* *ext) const 
    {
      double w[4];
      double result = 0.;

      if(!(*this->cacheReady))
      {
        for (int point_i = 0; point_i < n; point_i++) 
        {
          {
            w[0] = ext[0]->val[point_i];
            w[1] = ext[1]->val[point_i];
            w[2] = ext[2]->val[point_i];
            w[3] = ext[3]->val[point_i];

            double e_1_1[4] = {1, 0, 0, 0};
            double e_2_1[4] = {0, 1, 0, 0};
            double e_3_1[4] = {0, 0, 1, 0};
            double e_4_1[4] = {0, 0, 0, 1};

            num_flux->P_plus(this->P_plus_cache[point_i], w, e_1_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 4, w, e_2_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 8, w, e_3_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(this->P_plus_cache[point_i] + 12, w, e_4_1, e->nx[point_i], e->ny[point_i]);

            w[0] = ext[0]->val_neighbor[point_i];
            w[1] = ext[1]->val_neighbor[point_i];
            w[2] = ext[2]->val_neighbor[point_i];
            w[3] = ext[3]->val_neighbor[point_i];

            double e_1_2[4] = {1, 0, 0, 0};
            double e_2_2[4] = {0, 1, 0, 0};
            double e_3_2[4] = {0, 0, 1, 0};
            double e_4_2[4] = {0, 0, 0, 1};

            num_flux->P_minus(this->P_minus_cache[point_i], w, e_1_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 4, w, e_2_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 8, w, e_3_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(this->P_minus_cache[point_i] + 12, w, e_4_2, e->nx[point_i], e->ny[point_i]);
          }
        }
        *(const_cast<EulerInterfaceBilinearForm*>(this))->cacheReady = true;
      }

      int index = j * 4 + i;

      if(u->val == NULL)
        if(v->val == NULL)
          for (int point_i = 0; point_i < n; point_i++) 
            result -= wt[point_i] * (this->P_minus_cache[point_i][index] * u->val_neighbor[point_i]) * v->val_neighbor[point_i];
        else
          for (int point_i = 0; point_i < n; point_i++) 
            result += wt[point_i] * (this->P_minus_cache[point_i][index] * u->val_neighbor[point_i]) * v->val[point_i];
      else
        if(v->val == NULL)
          for (int point_i = 0; point_i < n; point_i++) 
            result -= wt[point_i] * (this->P_plus_cache[point_i][index] * u->val[point_i]) * v->val_neighbor[point_i];
        else
          for (int point_i = 0; point_i < n; point_i++) 
            result += wt[point_i] * (this->P_plus_cache[point_i][index] * u->val[point_i]) * v->val[point_i];
        
      return result * wf->get_current_time_step();
    }

    MatrixFormDG<double>* EulerInterfaceBilinearForm::clone()  const
    { 
      EulerInterfaceBilinearForm* form = new EulerInterfaceBilinearForm(this->i, this->j, this->num_flux->kappa, this->fluxes, this->cacheReady, this->P_plus_cache, this->P_minus_cache);
      form->wf = this->wf;
      return form;
    }

/////////////------------WEAKFORM-------------


	EulerInterface::EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa)), riemann_invariants(new RiemannInvariants(kappa))
	{

    P_plus_cache_DG = new double*[13];
    P_minus_cache_DG = new double*[13];

    for(int coordinate_i = 0; coordinate_i < 13; coordinate_i++)
		{
			P_plus_cache_DG[coordinate_i] = new double[16];
      P_minus_cache_DG[coordinate_i] = new double[16];
		}

		for(int form_i = 0; form_i < 4; form_i++)
			for(int form_j = 0; form_j < 4; form_j++)
			{
			EulerInterfaceBilinearForm* formDG = new EulerInterfaceBilinearForm(form_i, form_j, kappa, euler_fluxes, &this->cacheReadyDG, this->P_plus_cache_DG, this->P_minus_cache_DG);
				  add_matrix_form_DG(formDG);
			}

    
    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));


	};


	EulerInterface ::~EulerInterface ()
	{
		delete euler_fluxes;
		delete riemann_invariants;	
	};
	
	    WeakForm<double>* EulerInterface::clone() const
    {
      const_cast<EulerInterface*>(this)->warned_nonOverride = false;
      return new EulerInterface(*this);
    }


