#include "definitions.h"

CustomWeakForm::CustomWeakForm(double time_step,double theta, Solution<double>* sln_prev_time, std::string inlet, Mesh* mesh, bool all, bool DG) : WeakForm<double>(1), mesh(mesh)
{
 this->set_ext(sln_prev_time);
  if(all)
  {
    add_matrix_form(new CustomMatrixFormVol(0, 0, time_step,theta));
    add_vector_form(new CustomVectorFormVol(0,time_step,theta));
    }
   add_matrix_form_surf(new CustomMatrixFormSurface(0, 0));    
   if(DG) add_matrix_form_DG(new CustomMatrixFormInterface(0, 0));

}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i])/time_step
    -wt[i] *theta* (u->val[i] *(v->dx[i] * (0.5- e->y[i]) + v->dy[i] * (e->x[i]-0.5) ));

  return result;
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVol::clone() const
{
  return new CustomWeakForm::CustomMatrixFormVol(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  Func<Scalar>* u_prev_time = ext[0];
  for (int i = 0; i < n; i++)
	result += wt[i] *  u_prev_time->val[i] * v->val[i]/time_step + (1.-theta)*wt[i] * ( u_prev_time->val[i] * (v->dx[i] * (0.5- e->y[i]) + v->dy[i] *  (e->x[i]-0.5)));

  return result;
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVol::clone() const
{
  return new CustomWeakForm::CustomVectorFormVol(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormSurface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                      Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
  {
 	Real v_x = (0.5- e->y[i]);
 	Real v_y = (e->x[i]-0.5); 
   Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
   result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->val[i], Scalar(0), a_dot_n) * v->val[i];
  }
  return result;
}

double CustomWeakForm::CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* CustomWeakForm::CustomMatrixFormSurface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormSurface(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormInterface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                        Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
 	 Real v_x = (0.5- e->y[i]);
 	 Real v_y = (e->x[i]-0.5); 
    Real a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
   result += wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v; 
      
  }
  return result;
}

double CustomWeakForm::CustomMatrixFormInterface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                  Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormInterface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                             Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormDG<double>* CustomWeakForm::CustomMatrixFormInterface::clone() const
{
  return new CustomWeakForm::CustomMatrixFormInterface(*this);
}

double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    double v_x = (0.5- e->y[i]); double v_y =(e->x[i]-0.5);
    double a_dot_n = static_cast<CustomWeakForm*>(wf)->calculate_a_dot_v(v_x, v_y, e->nx[i], e->ny[i]);
    // Function values for Dirichlet boundary conditions.
    result += -wt[i] * static_cast<CustomWeakForm*>(wf)->upwind_flux(0, g<double,double>(static_cast<CustomWeakForm*>(wf)->mesh->get_boundary_markers_conversion().get_user_marker(e->edge_marker).marker), a_dot_n) * v->val[i];
  }
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += -wt[i] * v->val[i];
  return result;
}

VectorFormSurf<double>* CustomWeakForm::CustomVectorFormSurface::clone() const
{
  return new CustomWeakForm::CustomVectorFormSurface(*this);
}

template<typename Real>
Real CustomWeakForm::CustomVectorFormSurface::F(Real x, Real y) const
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormSurface::g(std::string ess_bdy_marker) const
{
  if(ess_bdy_marker == inlet) return 1; else return 0;
}

double CustomWeakForm::calculate_a_dot_v(double x, double y, double vx, double vy) const
{
 
 return  x*vx + y*vy;
}

Ord CustomWeakForm::calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const
{
  return Ord(10);
}

double CustomWeakForm::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomWeakForm::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}









//---------------Massematrix-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step,Solution<double>* sln_prev_time) : WeakForm<double>(1) {
 this->set_ext(sln_prev_time);
		CustomMatrixFormVolMassmatrix* mass_form= new CustomMatrixFormVolMassmatrix(0, 0, time_step);	
		add_matrix_form(mass_form);
		 

  }

    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolMassmatrix::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
		     Scalar result = Scalar(0); 
	  for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->val[i])/time_step;
	  return result;

    };

   double CustomMatrixFormVolMassmatrix::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

    Ord CustomMatrixFormVolMassmatrix::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

MatrixFormVol<double>* CustomMatrixFormVolMassmatrix::clone() const
{
  return new CustomMatrixFormVolMassmatrix(*this);
}



//---------------Konvektion-----------

  CustomWeakFormConvection::CustomWeakFormConvection(Solution<double>* sln_prev_time) : WeakForm<double>(1) {
   this->set_ext(sln_prev_time);
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
   VectorFormVolConvection* vector_form = new VectorFormVolConvection(0);
   add_vector_form(vector_form);
  };



    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
  result += wt[i] * (u->val[i] *(v->dx[i] * (0.5- e->y[i]) + v->dy[i] * (e->x[i]-0.5) ));
    //result += -wt[i] * (v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));
  return result;

    };

double CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* CustomMatrixFormVolConvection::clone() const
{
  return new CustomMatrixFormVolConvection(*this);
}




    template<typename Real, typename Scalar>
    Scalar VectorFormVolConvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
  Scalar result = Scalar(0); 
  Func<Scalar>* u_prev_time = ext[0];
  for (int i = 0; i < n; i++)
    result += -wt[i] * ( v->val[i] * (u_prev_time->dx[i] * (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
  return result;
    };
     double VectorFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

     Ord VectorFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };

   VectorFormVol<double>* VectorFormVolConvection::clone() const{
 			 return new VectorFormVolConvection(*this);

		}


//---------Gradient Reconstruction---------------
//R_H^1
 GradientReconstruction_1::GradientReconstruction_1( Solution<double>* sln) : WeakForm<double>(1) {
  this->set_ext(sln);
    add_matrix_form(new GradientReconstructionMatForm_1(0, 0));
    GradientReconstructionVectorForm_1* vector_form = new GradientReconstructionVectorForm_1(0);
    add_vector_form(vector_form);
  };


    template<typename Real, typename Scalar>
    Scalar GradientReconstructionMatForm_1 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

		       Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] );	
		return result;

    };

double GradientReconstructionMatForm_1 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };
Ord GradientReconstructionMatForm_1 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* GradientReconstructionMatForm_1::clone() const
{
  return new GradientReconstructionMatForm_1(*this);
}

    template<typename Real, typename Scalar>
    Scalar GradientReconstructionVectorForm_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
			Scalar result = Scalar(0);
			Func<Scalar>* u_h = ext[0];
			for (int i = 0; i < n; i++) 	
				result += wt[i]*(v->val[i]*u_h->dx[i]);	
			return result;

    };

  double GradientReconstructionVectorForm_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord GradientReconstructionVectorForm_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };
   VectorFormVol<double>* GradientReconstructionVectorForm_1::clone() const{
 			 return new GradientReconstructionVectorForm_1(*this);

		}
//R_H^2

 GradientReconstruction_2::GradientReconstruction_2( Solution<double>* sln) : WeakForm<double>(1) {
   this->set_ext(sln);
    add_matrix_form(new GradientReconstructionMatForm_2(0, 0));
    GradientReconstructionVectorForm_2* vector_form = new GradientReconstructionVectorForm_2(0);
    add_vector_form(vector_form);
  };



    template<typename Real, typename Scalar>
    Scalar GradientReconstructionMatForm_2 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

		       Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] );	
		return result;

    };

double GradientReconstructionMatForm_2 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };
Ord GradientReconstructionMatForm_2 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   MatrixFormVol<double>* GradientReconstructionMatForm_2::clone() const
{
  return new GradientReconstructionMatForm_2(*this);
}

    template<typename Real, typename Scalar>
    Scalar GradientReconstructionVectorForm_2::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
			Scalar result = Scalar(0);
			Func<Scalar>* u_h = ext[0];			
			for (int i = 0; i < n; i++) 				
				result += wt[i]*(v->val[i]*u_h->dy[i]);				
				
			return result;

    };

  double GradientReconstructionVectorForm_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord GradientReconstructionVectorForm_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };

   VectorFormVol<double>* GradientReconstructionVectorForm_2::clone() const{
 			 return new GradientReconstructionVectorForm_2(*this);

		}
//------------------- Initial condition ----------------

 void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const {
      
   	double radius = 0.;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) {		
		dx = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
		dy = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
	}
	/*else{			
		//cone
		x_0 = 0.5;
		y_0 = 0.25;
		radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
		if((radius< 1.0)&&(x!=x_0)) { 	
				dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
			dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
		}*/
		else{dx=0.; dy=0.;
		}	
  
	//}
		

};

 double CustomInitialCondition::value(double x, double y) const {
       
     double result = 0.0;
	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) { 
		 result = (1.0+ std::cos(PI*radius))/4.0;
		return result;	
	}
	//slotted cylinder
/*	x_0 = 0.5;
	y_0 = 0.75;
	radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if(radius <= 1) { 	
		if(fabs((x-x_0))>= 0.025) return 1.0;
		if(y>=0.85) return 1.0;
	}	
	//cone
	x_0 = 0.5;
	y_0 = 0.25;
	radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if(radius<= 1.0) { 	
		result = 1.0-radius;
	}	*/
       return result;
};

 Ord CustomInitialCondition::ord(Ord x, Ord y) const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition::clone() const
    {
return new CustomInitialCondition(this->mesh);

    }



