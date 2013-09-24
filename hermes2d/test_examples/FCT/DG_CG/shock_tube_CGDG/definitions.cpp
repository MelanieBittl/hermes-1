#include "definitions.h"


  EulerEquationsWeakForm_Mass::EulerEquationsWeakForm_Mass(double time_step, int num_of_equations): WeakForm<double>(num_of_equations), time_step(time_step), num_of_equations(num_of_equations)
	{
    add_matrix_form(new EulerEquationsBilinearFormTime(0));  //density
    add_matrix_form(new EulerEquationsBilinearFormTime(1));		//density_vel_x
    add_matrix_form(new EulerEquationsBilinearFormTime(2));		//density_vel_y
    add_matrix_form(new EulerEquationsBilinearFormTime(3));		//energy
	}


/*	EulerEquationsWeakForm_Mass::~EulerEquationsWeakForm_Mass()
		{
			for(int i=0; i<this->mfvol.size();i++){
				delete get_mfvol()[i];
			}
			WeakForm<double>::delete_all();
		}*/

	    WeakForm<double>* EulerEquationsWeakForm_Mass::clone() const
    {
      const_cast<EulerEquationsWeakForm_Mass*>(this)->warned_nonOverride = false;
      return new EulerEquationsWeakForm_Mass(this->time_step, this->num_of_equations);
    }
    
    
    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormTime::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double EulerEquationsBilinearFormTime::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormTime::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormTime::clone() const { return new EulerEquationsBilinearFormTime(component_i); }






//----------Matrix K--------------
  EulerEquationsWeakForm_K::EulerEquationsWeakForm_K(double kappa,double time_step, MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x, MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, int num_of_equations): WeakForm<double>(num_of_equations), euler_fluxes(new EulerFluxes(kappa))
	{
    add_matrix_form(new EulerEquationsBilinearFormDensity(0));
    add_matrix_form(new EulerEquationsBilinearFormDensity(1));
    add_matrix_form(new EulerEquationsBilinearFormDensity(2));
    add_matrix_form(new EulerEquationsBilinearFormDensity(3));

   add_matrix_form(new EulerEquationsBilinearFormDensityVelX(kappa,0));
   add_matrix_form(new EulerEquationsBilinearFormDensityVelX(kappa,1));
   add_matrix_form(new EulerEquationsBilinearFormDensityVelX(kappa,2));
   add_matrix_form(new EulerEquationsBilinearFormDensityVelX(kappa,3));

    add_matrix_form(new EulerEquationsBilinearFormDensityVelY(kappa,0));
    add_matrix_form(new EulerEquationsBilinearFormDensityVelY(kappa,1));
    add_matrix_form(new EulerEquationsBilinearFormDensityVelY(kappa,2));
    add_matrix_form(new EulerEquationsBilinearFormDensityVelY(kappa,3));

    add_matrix_form(new EulerEquationsBilinearFormEnergy(kappa,0));
    add_matrix_form(new EulerEquationsBilinearFormEnergy(kappa,1));
    add_matrix_form(new EulerEquationsBilinearFormEnergy(kappa,2));
    add_matrix_form(new EulerEquationsBilinearFormEnergy(kappa,3));

    //for(unsigned int vector_form_i = 0;vector_form_i < this->mfvol.size();vector_form_i++)     
     //mfvol.at(vector_form_i)->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    this->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
	};


	
/*	EulerEquationsWeakForm_K::~EulerEquationsWeakForm_K()
	{		
		for(int i=0; i<this->mfvol.size();i++){
			delete get_mfvol()[i];			
		}
		delete euler_fluxes;

		WeakForm<double>::delete_all();

	};*/
	    WeakForm<double>* EulerEquationsWeakForm_K::clone() const
    {
     // const_cast<EulerEquationsWeakForm_K*>(this)->warned_nonOverride = false;
    //  return new EulerEquationsWeakForm_K(*this);


   EulerEquationsWeakForm_K* wf;
    wf = new EulerEquationsWeakForm_K(*this);

    wf->ext.clear();

    for(unsigned int i = 0; i < this->ext.size(); i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->ext[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        wf->ext.push_back(new Solution<double>());
        wf->ext.back()->copy(this->ext[i]);
      }
      else
        wf->ext.push_back(this->ext[i]->clone());
    }
    return wf;



    }

    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensity::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
		 {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==1){
				
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if(j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }
    
     

    double EulerEquationsBilinearFormDensity::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    } 


    Ord EulerEquationsBilinearFormDensity::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensity::clone() const
    {
    return new EulerEquationsBilinearFormDensity(this->j);
    }



    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensityVelX::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const  
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerEquationsBilinearFormDensityVelX::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormDensityVelX::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensityVelX::clone() const { return new EulerEquationsBilinearFormDensityVelX(kappa,this->j); }



    template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormDensityVelY::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			if(j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerEquationsBilinearFormDensityVelY::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext,u, v, e, ext);
    }

    Ord EulerEquationsBilinearFormDensityVelY::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return Ord(20);
    }

    MatrixFormVol<double>* EulerEquationsBilinearFormDensityVelY::clone() const { return new EulerEquationsBilinearFormDensityVelY(kappa,this->j); }


  template<typename Real, typename Scalar>
    Scalar EulerEquationsBilinearFormEnergy::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
			 if (j==0){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dy[i];
			}else if (j==1){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
			}else if (j==2){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], ext[3]->val[i]) 
          * v->dy[i];
			}else if (j==3){
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * u->val[i] 
        * (static_cast<EulerEquationsWeakForm_K*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext[0]->val[i], ext[1]->val[i], ext[2]->val[i], Scalar(0)) 
          * v->dy[i];
				}
      }
      return result;
    }

    double EulerEquationsBilinearFormEnergy::value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u,v, e, ext);
    }

    Ord EulerEquationsBilinearFormEnergy::ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const
    {
      return Ord(20);
    }


    MatrixFormVol<double>* EulerEquationsBilinearFormEnergy::clone() const { return new EulerEquationsBilinearFormEnergy(kappa,this->j); }




