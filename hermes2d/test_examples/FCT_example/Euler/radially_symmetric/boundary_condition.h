#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include "euler_util.h"
#include "hermes2d.h"
 #define PI (3.141592653589793)   


//Boundary condition are imposed in the weak sense


//Boundary Flux Jacobians - > Matrix dS
 class EulerBoundaryBilinearForm_rho: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_rho(double kappa, int j) : MatrixFormSurf<double>(0,j), kappa(kappa),j(j) {}

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

 class EulerBoundaryBilinearForm_vel_x: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_x(double kappa,int j) : MatrixFormSurf<double>(1,j), kappa(kappa),j(j) {}

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

 class EulerBoundaryBilinearForm_vel_y: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_y(double kappa,int j) : MatrixFormSurf<double>(2,j), kappa(kappa),j(j) {}

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
 class EulerBoundaryBilinearForm_e: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_e(double kappa,int j) : MatrixFormSurf<double>(3,j), kappa(kappa),j(j) {}

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

//--------Linearform for Surf---------

 class EulerBoundary_rho: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_rho(double kappa) : VectorFormSurf<double>(0), kappa(kappa) {}

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

 class EulerBoundary_v_x: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_v_x(double kappa) : VectorFormSurf<double>(1), kappa(kappa) {}

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
 class EulerBoundary_v_y: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_v_y(double kappa) : VectorFormSurf<double>(2), kappa(kappa) {}

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

 class EulerBoundary_e: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_e(double kappa) : VectorFormSurf<double>(3), kappa(kappa) {}

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


//WeakForm for boundary ->Matrix dS

class EulerBoundary : public WeakForm<double>
{
public:

  EulerBoundary(double kappa,Solution<double>* rho_ext, Solution<double>* v1_ext, Solution<double>* v2_ext, Solution<double>* energy_ext, 
Solution<double>* prev_density, Solution<double>* prev_density_vel_x,  Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, 
bool mirror_condition= true ,
int num_of_equations = 4);

	~EulerBoundary();

  // Members.
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;


};









#endif
