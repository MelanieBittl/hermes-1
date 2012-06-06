#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "euler_util.h"
#include "hermes2d.h"
 #define PI (3.141592653589793)   



using namespace Hermes;
using namespace Hermes::Hermes2D;



//-----------------------------Bilinearform for Time-Discretization

 class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const ;

   virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const ;

    MatrixFormVol<double>* clone();
  };


  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i), component_i(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const ;


    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const ;

    VectorFormVol<double>* clone() ;
    // Member.
    int component_i;
  };
//----------------Bilinearform for Surf

 class EulerEquationsBilinearFormBoundary_rho: public MatrixFormSurf<double>
  {
  public:
    EulerEquationsBilinearFormBoundary_rho(double kappa, int j) : MatrixFormSurf<double>(0,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormSurf<double>* clone();

		double kappa;
	int j;
	};

 class EulerEquationsBilinearFormBoundary_vel_x: public MatrixFormSurf<double>
  {
  public:
    EulerEquationsBilinearFormBoundary_vel_x(double kappa,int j) : MatrixFormSurf<double>(1,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormSurf<double>* clone();

		double kappa;
		int j;
	};

 class EulerEquationsBilinearFormBoundary_vel_y: public MatrixFormSurf<double>
  {
  public:
    EulerEquationsBilinearFormBoundary_vel_y(double kappa,int j) : MatrixFormSurf<double>(2,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormSurf<double>* clone();

		double kappa;
		int j;
	};
 class EulerEquationsBilinearFormBoundary_e: public MatrixFormSurf<double>
  {
  public:
    EulerEquationsBilinearFormBoundary_e(double kappa,int j) : MatrixFormSurf<double>(3,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormSurf<double>* clone();

		double kappa;
		int j;
	};

//--------Linearform---------

 class EulerEquationsLinearFormBoundary_rho: public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormBoundary_rho(double kappa) : VectorFormSurf<double>(0), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    VectorFormSurf<double>* clone();

		double kappa;
	};

 class EulerEquationsLinearFormBoundary_v_x: public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormBoundary_v_x(double kappa) : VectorFormSurf<double>(1), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    VectorFormSurf<double>* clone();

		double kappa;
	};
 class EulerEquationsLinearFormBoundary_v_y: public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormBoundary_v_y(double kappa) : VectorFormSurf<double>(2), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    VectorFormSurf<double>* clone();

		double kappa;
	};

 class EulerEquationsLinearFormBoundary_e: public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormBoundary_e(double kappa) : VectorFormSurf<double>(3), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    VectorFormSurf<double>* clone();

		double kappa;
	};

//----------------Bilinearform for K

 class EulerEquationsBilinearFormDensity : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensity(int j) : MatrixFormVol<double>(0,j), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormVol<double>* clone();
		int j;
	};

 class EulerEquationsBilinearFormDensityVelX : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensityVelX(double kappa, int j)
      : MatrixFormVol<double>(1,j), kappa(kappa), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;


    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const ;


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const ;

   MatrixFormVol<double>* clone();

    double kappa;
		int j;
  };

  class EulerEquationsBilinearFormDensityVelY : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensityVelY(double kappa,int j) 
      : MatrixFormVol<double>(2,j), kappa(kappa), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;
   

    double value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const; 

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const; 

    MatrixFormVol<double>* clone();

    double kappa;
		int j;
  };


 class EulerEquationsBilinearFormEnergy : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormEnergy(double kappa, int j) 
      : MatrixFormVol<double>(3,j), kappa(kappa), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[],Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const ;


    double value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const ;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const;

   MatrixFormVol<double>* clone();

    double kappa;
		int j;
  };




//---------------

class EulerEquationsWeakForm_Mass : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_Mass(double time_step, Solution<double>* prev_density, Solution<double>* prev_density_vel_x, 
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations = 4);
	~EulerEquationsWeakForm_Mass();

};

class EulerEquationsWeakForm_K : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_K(double kappa,double time_step, Solution<double>* prev_density, Solution<double>* prev_density_vel_x, 
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations = 4);
	~EulerEquationsWeakForm_K();

  // Members.
  EulerFluxes* euler_fluxes;

};


class EulerEquationsWeakForm_Surf : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_Surf(double kappa,Solution<double>* prev_density, Solution<double>* prev_density_vel_x, 
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, int num_of_equations = 4);
	~EulerEquationsWeakForm_Surf();

  // Members.
  EulerFluxes* euler_fluxes;


};



#endif

