#include "definitions.h"
//---------------Massematrix-----------
 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step,Solution<double>* sln_prev_time) : WeakForm<double>(1) {
		CustomMatrixFormVolMassmatrix* mass_form= new CustomMatrixFormVolMassmatrix(0, 0, time_step);	
		add_matrix_form(mass_form);
		VectorFormVolMass* vector_form = new VectorFormVolMass(0, time_step);
		vector_form->ext.push_back(sln_prev_time);
		add_vector_form(vector_form);
  }
 CustomWeakFormMassmatrix::~CustomWeakFormMassmatrix(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();

	};


    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolMassmatrix::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
		     Scalar result = Scalar(0); 
	  for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->val[i])/time_step;
	  return result;

    };

   double CustomMatrixFormVolMassmatrix::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

    Ord CustomMatrixFormVolMassmatrix::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };



    template<typename Real, typename Scalar>
    Scalar VectorFormVolMass::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
	  Scalar result = Scalar(0);
	  Func<Scalar>* u_prev_time = ext->fn[0];
	  for (int i = 0; i < n; i++)
		result += wt[i] *  u_prev_time->val[i] * v->val[i]/time_step;
	  return result;

    };

   double VectorFormVolMass::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord VectorFormVolMass::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };

//---------------Konvektion-----------

  CustomWeakFormConvection::CustomWeakFormConvection(Solution<double>* sln_prev_time) : WeakForm<double>(1) {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
   VectorFormVolConvection* vector_form = new VectorFormVolConvection(0);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };
	CustomWeakFormConvection::~CustomWeakFormConvection(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();
	};




    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

     Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++)
    result += -wt[i] * (v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));
  return result;

    };

double CustomMatrixFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   



    template<typename Real, typename Scalar>
    Scalar VectorFormVolConvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = Scalar(0); 
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += -wt[i] * ( v->val[i] * (u_prev_time->dx[i] * (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
  return result;
    };
     double VectorFormVolConvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

     Ord VectorFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };


//------Matrix & Vektorform for higher Order solution----------------


  ConvectionForm::ConvectionForm( double tau, Solution<double>* sln_prev_time) : WeakForm<double>(1) {
    add_matrix_form(new ConvectionMatForm(0, 0,  tau));
    VectorConvection* vector_form = new VectorConvection(0,  tau);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  }
 ConvectionForm::~ConvectionForm(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();
	};


 template<typename Real, typename Scalar>
    Scalar ConvectionMatForm::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
         Scalar result = Scalar(0); 
				for (int i = 0; i < n; i++) 
					result += wt[i] * (u->val[i] * v->val[i] / tau + 0.5*v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));	
				return result;

    }

   double ConvectionMatForm::value(int n, double *wt, Func<double>  *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double>  *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

   Ord ConvectionMatForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    template<typename Real, typename Scalar>
    Scalar VectorConvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = Scalar(0);
  Func<Scalar>* u_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) 
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau +
                      0.5*v->val[i] * (u_prev_newton->dx[i] * (0.5- e->y[i]) + u_prev_newton->dy[i] *  (e->x[i]-0.5) +
			u_prev_time->dx[i] *  (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
	
  return result;

    }

    double VectorConvection::value(int n, double *wt, Func<double > *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

     Ord VectorConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }



//-----------------Residual for error estimator--------------------------------

 ResidualForm::ResidualForm( double tau, Solution<double>* sln_prev_time, Solution<double>* ref) : WeakForm<double>(1) {
    add_matrix_form(new ResidualMatForm(0, 0,  tau));
    VectorResidual* vector_form = new VectorResidual(0,  tau);
    vector_form->ext.push_back(sln_prev_time);
		vector_form->ext.push_back(ref);
    add_vector_form(vector_form);
  };

	ResidualForm::~ResidualForm(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();
	};

    template<typename Real, typename Scalar>
    Scalar ResidualMatForm::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

		       Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] / tau + 0.5*v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));	
		return result;

    };

double ResidualMatForm::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };
Ord ResidualMatForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };



    template<typename Real, typename Scalar>
    Scalar VectorResidual::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
			Scalar result = Scalar(0);
			Func<Scalar>* u_prev_time = ext->fn[0];
			Func<Scalar>* u_h = ext->fn[1];
			for (int i = 0; i < n; i++) 	
				result += wt[i]*(u_prev_time->val[i] * v->val[i] / tau -
						0.5*v->val[i]*(	u_prev_time->dx[i] *  (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5))-						
								u_h->val[i] * v->val[i] / tau - 
							(0.5*v->val[i]*(u_h->dx[i] * (0.5- e->y[i]) + u_h->dy[i] * (e->x[i]-0.5))));							
	
			return result;

    };

  double VectorResidual::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord VectorResidual::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };


//---------Gradient Reconstruction---------------
//R_H^1
 GradientReconstruction_1::GradientReconstruction_1( Solution<double>* sln) : WeakForm<double>(1) {
    add_matrix_form(new GradientReconstructionMatForm_1(0, 0));
    GradientReconstructionVectorForm_1* vector_form = new GradientReconstructionVectorForm_1(0);
    vector_form->ext.push_back(sln);
    add_vector_form(vector_form);
  };

	GradientReconstruction_1::~GradientReconstruction_1(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();
	};

    template<typename Real, typename Scalar>
    Scalar GradientReconstructionMatForm_1 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

		       Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] );	
		return result;

    };

double GradientReconstructionMatForm_1 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };
Ord GradientReconstructionMatForm_1 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };



    template<typename Real, typename Scalar>
    Scalar GradientReconstructionVectorForm_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
			Scalar result = Scalar(0);
			Func<Scalar>* u_h = ext->fn[0];
			for (int i = 0; i < n; i++) 	
				result += wt[i]*(v->val[i]*u_h->dx[i]);	
			return result;

    };

  double GradientReconstructionVectorForm_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord GradientReconstructionVectorForm_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };

//R_H^2

 GradientReconstruction_2::GradientReconstruction_2( Solution<double>* sln) : WeakForm<double>(1) {
    add_matrix_form(new GradientReconstructionMatForm_2(0, 0));
    GradientReconstructionVectorForm_2* vector_form = new GradientReconstructionVectorForm_2(0);
    vector_form->ext.push_back(sln);
    add_vector_form(vector_form);
  };

	GradientReconstruction_2::~GradientReconstruction_2(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm<double>::delete_all();
	};

    template<typename Real, typename Scalar>
    Scalar GradientReconstructionMatForm_2 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

		       Scalar result = Scalar(0); 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] );	
		return result;

    };

double GradientReconstructionMatForm_2 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    };
Ord GradientReconstructionMatForm_2 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };



    template<typename Real, typename Scalar>
    Scalar GradientReconstructionVectorForm_2::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
			Scalar result = Scalar(0);
			Func<Scalar>* u_h = ext->fn[0];
			for (int i = 0; i < n; i++) 	
				result += wt[i]*(v->val[i]*u_h->dy[i]);	
			return result;

    };

  double GradientReconstructionVectorForm_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    };

    Ord GradientReconstructionVectorForm_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };


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
	else{			
		//cone
		x_0 = 0.5;
		y_0 = 0.25;
		radius = 1.0/0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0));
		if((radius< 1.0)&&(x!=x_0)) { 	
				dx = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*x;
			dy = -(1.0/(0.15 * std::sqrt( std::pow((x-x_0),2.0) + std::pow((y-y_0),2.0))))*2*y;	
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

 Ord CustomInitialCondition::ord(Ord x, Ord y) const {
      return Ord(10);
};





