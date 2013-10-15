#include "definitions.h"


double calc_abs_v(Element* e)
{
Hermes::Hermes2D::ElementMode2D mode = HERMES_MODE_QUAD;
			if(e->is_triangle()) mode = HERMES_MODE_TRIANGLE;
	// order of integral
	const int order = 4;
	// refmap for computing Jacobian
	RefMap* rm = new RefMap;
	rm->set_quad_2d(&g_quad_2d_std);

			// get the quadrature points
			int np = g_quad_2d_std.get_num_points(order,mode);
			double3 *pt = g_quad_2d_std.get_points(order, mode);
			// get the constant Jacobian
		  rm->set_active_element(e);

		//double diam = e->get_diameter();
		//double area =Hermes::sqrt(e->get_area());

		double abs_v =0;	
		double* x_coord = rm->get_phys_x(order);
		double* y_coord = rm->get_phys_y(order);
			for( int j = 0; j < np; ++j )
			{
				double v_x = 0.5-y_coord[j];
				double v_y = x_coord[j]-0.5;	 
				abs_v += pt[j][2]*(v_x*v_x+v_y*v_y);
			}
	
	delete rm;
	return Hermes::sqrt(abs_v);
}

//---------------Massematrix-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step,MeshFunctionSharedPtr<double> sln_prev_time) : WeakForm<double>(1) {
 this->set_ext(sln_prev_time);
		CustomMatrixFormVolMassmatrix* mass_form= new CustomMatrixFormVolMassmatrix(0, 0, time_step);	
		add_matrix_form(mass_form);
		VectorFormVolMass* vector_form = new VectorFormVolMass(0, time_step);		
		add_vector_form(vector_form);
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

    template<typename Real, typename Scalar>
    Scalar VectorFormVolMass::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
	  Scalar result = Scalar(0);
	  Func<Scalar>* u_prev_time = ext[0];
	  for (int i = 0; i < n; i++)
		result += wt[i] *  u_prev_time->val[i] * v->val[i]/time_step;
	  return result;

    };

   double VectorFormVolMass::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord VectorFormVolMass::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };
   VectorFormVol<double>* VectorFormVolMass::clone() const{
 			 return new VectorFormVolMass(*this);

		}
//---------------Konvektion-----------

//----------------------wf-

  CustomWeakFormConvection::CustomWeakFormConvection() : WeakForm<double>(1) {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
  };


CustomWeakFormInterface::CustomWeakFormInterface(): WeakForm<double>(1) {
	 add_matrix_form_DG(new CustomMatrixFormInterface(0, 0)); 
  };


CustomWeakFormSurface::CustomWeakFormSurface(): WeakForm<double>(1) {
	 add_matrix_form_surf(new CustomMatrixFormSurface(0, 0)); 
  };

//------------surface---------
double CustomMatrixFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                                Geom<double> *e, Func<double> **ext) const
{
  double result = 0.;
	for (int i = 0; i < n; i++)
	{
			double v_x =(0.5- e->y[i]);
			double v_y = (e->x[i]-0.5) ; 
   double a_dot_n = (v_x*e->nx[i]+v_y*e->ny[i]);
	if(a_dot_n>=0)
   	result -= wt[i] *u->val[i]*a_dot_n * v->val[i];
		
	}
		return result;

}

Ord CustomMatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                           Geom<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val[i]*u->val[i];
  return result;
}

MatrixFormSurf<double>* CustomMatrixFormSurface::clone() const
{
  return new CustomMatrixFormSurface(*this);
}


//---------------convection------------
    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
{		Real v_x = (0.5- e->y[i]); 
 		Real v_y = (e->x[i]-0.5) ; 
		result += wt[i]*u->val[i]*(v->dx[i]*v_x+ v->dy[i]*v_y);

}
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

//-----------interface-------------------------
template<typename Real, typename Scalar>
Scalar CustomMatrixFormInterface::matrix_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, DiscontinuousFunc<Real> *u, DiscontinuousFunc<Real> *v, Geom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
		Real v_x = (0.5- e->y[i]); 
 		Real v_y = (e->x[i]-0.5) ; 
    Real a_dot_n = (v_x*e->nx[i]+ v_y*e->ny[i]);
    Real jump_v = (v->fn_central == NULL ? -v->val_neighbor[i] : v->val[i]);
    if(u->fn_central == NULL)
      result += wt[i] * upwind_flux(Real(0), u->val_neighbor[i], a_dot_n) * jump_v;
    else
      result += wt[i] * upwind_flux(u->val[i], Real(0), a_dot_n) * jump_v;
      
  }
  return (-result);
}

double CustomMatrixFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomMatrixFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormDG<double>* CustomMatrixFormInterface::clone() const
{
  return new CustomMatrixFormInterface(*this);
}

double CustomMatrixFormInterface::upwind_flux(double u_cent, double u_neib, double a_dot_n) const
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib);
}

Ord CustomMatrixFormInterface::upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const
{
  return a_dot_n * (u_cent + u_neib);
}

//---------Streamline Diffusion---------------

CustomWeakFormStreamlineDiffusion::CustomWeakFormStreamlineDiffusion (MeshSharedPtr mesh): WeakForm<double>(1) {
	 add_matrix_form(new Streamline(0,0,mesh)); 
  };

CustomWeakFormStreamlineDiffusionMass::CustomWeakFormStreamlineDiffusionMass (MeshSharedPtr mesh, double ts): WeakForm<double>(1) {
	 add_matrix_form(new StreamlineMass(0,0,mesh,ts)); 
  };


template<typename Real, typename Scalar>
Scalar Streamline::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);

		Real diam = e->diam;
		//Real area = Hermes::sqrt(e->area);
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem);
		for (int i = 0; i < n; i++)
		{
				Real v_x = (0.5- e->y[i]); 
				Real v_y = (e->x[i]-0.5) ; 
				Real tau = diam/(2.*abs_v);
				result -= wt[i] *(tau*(u->dx[i] * v_x + u->dy[i] * v_y ) *(v->dx[i] * v_x + v->dy[i] * v_y ));

		}


  return result;
}

double Streamline::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord Streamline::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* Streamline::clone() const
{
  return new Streamline(*this);
}


template<typename Real, typename Scalar>
Scalar StreamlineMass::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                                                  Geom<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);

		Real diam = e->diam;
		//Real area = Hermes::sqrt(e->area);
		Element* elem = mesh->get_element(e->id);
		double abs_v =calc_abs_v(elem);
		for (int i = 0; i < n; i++)
		{
				Real v_x = (0.5- e->y[i]); 
				Real v_y = (e->x[i]-0.5) ; 
				Real tau = diam/(2.*abs_v);
				result += wt[i] *(tau*(u->val[i]/ts) *(v->dx[i] * v_x + v->dy[i] * v_y ));
		}


  return result;
}

double StreamlineMass::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                                            Geom<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord StreamlineMass::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                       Geom<Ord> *e, Func<Ord> **ext) const
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* StreamlineMass::clone() const
{
  return new StreamlineMass(*this);
}


//---------Gradient Reconstruction---------------
//R_H^1
 GradientReconstruction_1::GradientReconstruction_1( MeshFunctionSharedPtr<double> sln) : WeakForm<double>(1) {
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

 GradientReconstruction_2::GradientReconstruction_2( MeshFunctionSharedPtr<double> sln) : WeakForm<double>(1) {
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
//double x_0 =0.5;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) {		
		dx = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(x-x_0);
		dy = -std::sin(radius*PI)/4.0*(PI/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(y-y_0);	
	}
	else{			
		//cone
		x_0 = 0.5;
		y_0 = 0.25;
		radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
		if((radius< 1.0)&&(x!=x_0)) { 	
				dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(x-x_0);
			dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2.*(x-x_0);	
		}
		else{dx=0.; dy=0.;
		}	
  
	}
		

};

 double CustomInitialCondition::value(double x, double y) const {
       
     double result = 0.0;
	double radius;
        //hump
	double x_0 =0.25;
//double x_0 =0.5;
	double y_0= 0.5;	
	radius = (1.0/0.15) * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
	if( radius<= 1.0) { 
		 result = (1.0+ std::cos(PI*radius))/4.0;
		return result;	
	}
	//slotted cylinder
x_0 = 0.5;
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
	}	
       return result;
};

 Ord CustomInitialCondition::ord(double x, double y)  const {
      return Ord(10);
};

 MeshFunction<double>* CustomInitialCondition::clone() const
    {
return new CustomInitialCondition(this->mesh);

    }



