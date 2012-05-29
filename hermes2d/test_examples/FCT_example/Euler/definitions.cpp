#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(double time_step, Solution<double>* prev_density, Solution<double>* prev_density_vel_x, 
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations): WeakForm<double>(num_of_equations)
{

    add_matrix_form(new EulerEquationsBilinearFormTime(0));  //density
    add_matrix_form(new EulerEquationsBilinearFormTime(1));		//density_vel_x
    add_matrix_form(new EulerEquationsBilinearFormTime(2));		//density_vel_y
    add_matrix_form(new EulerEquationsBilinearFormTime(3));		//energy

    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol.size();vector_form_i++) 
    {
      vfvol.at(vector_form_i)->ext.push_back(prev_density);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfvol.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfvol.at(vector_form_i)->ext.push_back(prev_energy);
    }
}


	EulerEquationsWeakForm_Mass::~EulerEquationsWeakForm_Mass()
	{
		for(int i=0; i<this->vfvol.size();i++){
			delete get_mfvol()[i];
			delete get_vfvol()[i];
		}
		WeakForm<double>::delete_all();


	}


    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormTime::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double EulerEquationsBilinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormTime::clone() { return new EulerEquationsBilinearFormTime(this->i); }


    template<typename Real, typename Scalar>
    Scalar EulerEquationsLinearFormTime::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, ext->fn[component_i], v);
    }

    double EulerEquationsLinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord EulerEquationsLinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    VectorFormVol<double>* EulerEquationsLinearFormTime::clone() { return new EulerEquationsLinearFormTime(*this); }

//-----------Surface--------------------------
  EulerEquationsWeakForm_Surf::EulerEquationsWeakForm_Surf(double kappa, Solution<double>* prev_density, Solution<double>* prev_density_vel_x, 
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa))
	{

 		add_matrix_form_surf(new EulerEquationsBilinearFormBoundary_rho(kappa));
    add_matrix_form_surf(new EulerEquationsBilinearFormBoundary_vel_x( kappa));
    add_matrix_form_surf(new EulerEquationsBilinearFormBoundary_vel_y( kappa));
    add_matrix_form_surf(new EulerEquationsBilinearFormBoundary_e( kappa));

  /*add_vector_form_surf(new EulerEquationsLinearFormBoundary_rho(kappa));
    add_vector_form_surf(new EulerEquationsLinearFormBoundary_vel_x( kappa));
    add_vector_form_surf(new EulerEquationsLinearFormBoundary_vel_y( kappa));
    add_vector_form_surf(new EulerEquationsLinearFormBoundary_e(kappa));*/

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf.size();vector_form_i++) 
    {
      vfsurf.at(vector_form_i)->ext.push_back(prev_density);
      vfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      vfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      vfsurf.at(vector_form_i)->ext.push_back(prev_energy);
    }
    for(unsigned int vector_form_i = 0;vector_form_i < this->mfsurf.size();vector_form_i++) 
    {
      mfsurf.at(vector_form_i)->ext.push_back(prev_density);
      mfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      mfsurf.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      mfsurf.at(vector_form_i)->ext.push_back(prev_energy);
    }

	};


	EulerEquationsWeakForm_Surf::~EulerEquationsWeakForm_Surf()
	{
		for(int i=0; i<this->vfsurf.size();i++)
		{
			delete get_mfsurf()[i];
			delete get_vfsurf()[i];
			}
		WeakForm<double>::delete_all();
		delete euler_fluxes;
	};


  template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormBoundary_rho::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const  
		 {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {		/*double n_x = 0.;
					double n_y = 0.;
					if(e->y[i]==-0.5) n_y = -1.;
					else if(e->y[i]==0.5) n_y = 1.;
					else if(e->x[i]==-0.5) n_x = -1.;
					else if(e->x[i]==0.5) n_x = 1.;*/

        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
      }
      return (-result);
    }  


    double EulerEquationsBilinearFormBoundary_rho::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerEquationsBilinearFormBoundary_rho::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerEquationsBilinearFormBoundary_rho::clone()
    {
					return new EulerEquationsBilinearFormBoundary_rho(*this);
    }


    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormBoundary_vel_x::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const  
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
      }
            return (-result);
    }

    double EulerEquationsBilinearFormBoundary_vel_x::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormBoundary_vel_x::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerEquationsBilinearFormBoundary_vel_x::clone() { return new EulerEquationsBilinearFormBoundary_vel_x(*this); }


  template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormBoundary_vel_y::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
      }
            return (-result);
    }

    double EulerEquationsBilinearFormBoundary_vel_y::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormBoundary_vel_y::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerEquationsBilinearFormBoundary_vel_y::clone() { return new EulerEquationsBilinearFormBoundary_vel_y(*this); }

  template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormBoundary_e::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
           * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * e->ny[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
           * e->nx[i]*v->val[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_Surf*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * e->ny[i]*v->val[i];
      }
            return (-result);
    }

    double EulerEquationsBilinearFormBoundary_e::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormBoundary_e::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* EulerEquationsBilinearFormBoundary_e::clone() { return new EulerEquationsBilinearFormBoundary_e(*this); }


//----------Matrix K--------------
  EulerEquationsWeakForm_K::EulerEquationsWeakForm_K(double kappa,double time_step, Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa))
	{
    add_matrix_form(new EulerEquationsBilinearFormDensity());
    add_matrix_form(new EulerEquationsBilinearFormDensityVelX(kappa));
    add_matrix_form(new EulerEquationsBilinearFormDensityVelY(kappa));
    add_matrix_form(new EulerEquationsBilinearFormEnergy(kappa));

    for(unsigned int vector_form_i = 0;vector_form_i < this->mfvol.size();vector_form_i++) 
    {
      mfvol.at(vector_form_i)->ext.push_back(prev_density);
      mfvol.at(vector_form_i)->ext.push_back(prev_density_vel_x);
      mfvol.at(vector_form_i)->ext.push_back(prev_density_vel_y);
      mfvol.at(vector_form_i)->ext.push_back(prev_energy);
    }
	};


	
	EulerEquationsWeakForm_K::~EulerEquationsWeakForm_K()
	{		
		for(int i=0; i<this->vfvol.size();i++){
			delete get_mfvol()[i];
			delete get_vfvol()[i];
		}
		delete euler_fluxes;

		WeakForm<double>::delete_all();

	};


    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensity::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const  
		 {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result;
    }
    
     

    double EulerEquationsBilinearFormDensity::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerEquationsBilinearFormDensity::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensity::clone()
    {
      EulerEquationsBilinearFormDensity* form = new EulerEquationsBilinearFormDensity(*this);
      form->wf = this->wf;
      return form;
    }



    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensityVelX::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const  
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result;
    }

    double EulerEquationsBilinearFormDensityVelX::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormDensityVelX::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensityVelX::clone() { return new EulerEquationsBilinearFormDensityVelX(*this); }



    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensityVelY::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result;
    }

    double EulerEquationsBilinearFormDensityVelY::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormDensityVelY::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensityVelY::clone() { return new EulerEquationsBilinearFormDensityVelY(*this); }


  template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormEnergy::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dy[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result;
    }

    double EulerEquationsBilinearFormEnergy::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormEnergy::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const
    {
      return Ord(20);
    }


    MatrixFormVol<double>* EulerEquationsBilinearFormEnergy::clone() { return new EulerEquationsBilinearFormEnergy(*this); }




