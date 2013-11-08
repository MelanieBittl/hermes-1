#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "euler_util.h"
#include "euler_flux.h"
#include "hermes2d.h"
 #define PI (3.141592653589793)   



using namespace Hermes;
using namespace Hermes::Hermes2D;



//-----------------------------Bilinearform for Time-Discretization

 class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i),component_i(i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const ;

   virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const ;

    MatrixFormVol<double>* clone() const;
    // Member.
    int component_i;
  };




//---------------WeakForms---------------
class EulerEquationsWeakForm_Mass : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_Mass(int num_of_equations = 4);

    WeakForm<double>* clone() const;

double time_step; int num_of_equations;
};



class EulerKS : public WeakForm<double>
{
public:

  EulerKS(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, 
MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, 
bool mirror_condition= true , int num_of_equations = 4);

	~EulerKS();

    WeakForm<double>* clone() const;


  // Members.
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	double kappa;

  MeshFunctionSharedPtr<double> prev_density;
  MeshFunctionSharedPtr<double> prev_density_vel_x;
  MeshFunctionSharedPtr<double> prev_density_vel_y;
  MeshFunctionSharedPtr<double> prev_energy;

  MeshFunctionSharedPtr<double> rho_ext;
  MeshFunctionSharedPtr<double> v1_ext;
  MeshFunctionSharedPtr<double> v2_ext;
  MeshFunctionSharedPtr<double> energy_ext;

protected:

 class EulerBoundaryBilinearForm_rho: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_rho(double kappa, int j) : MatrixFormSurf<double>(0,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

		double kappa;
	int j;
	};

 class EulerBoundaryBilinearForm_vel_x: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_x(double kappa,int j) : MatrixFormSurf<double>(1,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

		double kappa;
		int j;
	};

 class EulerBoundaryBilinearForm_vel_y: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_y(double kappa,int j) : MatrixFormSurf<double>(2,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

		double kappa;
		int j;
	};
 class EulerBoundaryBilinearForm_e: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_e(double kappa,int j) : MatrixFormSurf<double>(3,j), kappa(kappa),j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

		double kappa;
		int j;
	};

 class EulerEquationsBilinearFormDensity : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensity(int j) : MatrixFormVol<double>(0,j), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;
		int j;
	};

 class EulerEquationsBilinearFormDensityVelX : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensityVelX(double kappa, int j)
      : MatrixFormVol<double>(1,j), kappa(kappa), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;


    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const ;


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const ;

   MatrixFormVol<double>* clone() const;

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
      Geom<Real> *e, Func<Scalar>  **ext) const ;
   

    double value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;

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
      Geom<Real> *e, Func<Scalar>  **ext) const ;


    double value(int n, double *wt, Func<double> *u_ext[],Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const ;

    Ord ord(int n, double *wt, Func<Ord> *u_ext[],Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const;

   MatrixFormVol<double>* clone() const;

    double kappa;
		int j;
  };



};




#endif

