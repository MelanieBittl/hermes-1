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
    EulerBoundaryBilinearForm_rho(double kappa, int j,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : MatrixFormSurf<double>(0,j), kappa(kappa),j(j),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition) {}

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
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};

 class EulerBoundaryBilinearForm_vel_x: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_x(double kappa,int j,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : MatrixFormSurf<double>(1,j), kappa(kappa),j(j),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

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
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};

 class EulerBoundaryBilinearForm_vel_y: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_vel_y(double kappa,int j,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : MatrixFormSurf<double>(2,j), kappa(kappa),j(j),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

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
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};
 class EulerBoundaryBilinearForm_e: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm_e(double kappa,int j,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : MatrixFormSurf<double>(3,j), kappa(kappa),j(j),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

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
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};
//--------Linearform for Surf---------

 class EulerBoundary_rho: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_rho(double kappa,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : VectorFormSurf<double>(0), kappa(kappa),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double kappa;
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};

 class EulerBoundary_v_x: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_v_x(double kappa,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : VectorFormSurf<double>(1), kappa(kappa),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double kappa;
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};
 class EulerBoundary_v_y: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_v_y(double kappa,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : VectorFormSurf<double>(2), kappa(kappa),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double kappa;
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};

 class EulerBoundary_e: public VectorFormSurf<double>
  {
  public:
    EulerBoundary_e(double kappa,  EulerFluxes* euler_fluxes,	RiemannInvariants* riemann_invariants,	bool mirror_condition) : VectorFormSurf<double>(3), kappa(kappa),euler_fluxes(euler_fluxes), riemann_invariants(riemann_invariants), mirror_condition(mirror_condition)  {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double kappa;
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
	};


//WeakForm for boundary ->Matrix dS und vector S

class EulerBoundary : public WeakForm<double>
{
public:

  EulerBoundary(double kappa,Hermes::vector<MeshFunctionSharedPtr<double> > slns, bool mirror_condition= true , int num_of_equations = 4);

	~EulerBoundary();

    WeakForm<double>* clone() const;


  // Members.
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	bool mirror_condition;
double kappa;

Hermes::vector<MeshFunctionSharedPtr<double> > slns;


};









#endif
