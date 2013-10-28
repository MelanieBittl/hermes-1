#include "interface.h"

/////////////------------WEAKFORM-------------


	EulerInterface::EulerInterface(double kappa, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa))
	{

		for(int form_i = 0; form_i < 4; form_i++)
		{
				EulerEquationsVectorFormFlux* formDG = new EulerEquationsVectorFormFlux(form_i, kappa, euler_fluxes);
				  add_vector_form_DG(formDG);
		}

   this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

	};


	EulerInterface ::~EulerInterface ()
	{
		delete euler_fluxes;

	};
	
	    WeakForm<double>* EulerInterface::clone() const
    {
      const_cast<EulerInterface*>(this)->warned_nonOverride = false;
      return new EulerInterface(*this);
    }


