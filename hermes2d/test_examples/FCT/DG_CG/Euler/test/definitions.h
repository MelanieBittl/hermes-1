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



class EulerK : public WeakForm<double>
{
public:

  EulerK(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, 
MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, 
bool mirror_condition= true , int num_of_equations = 4);

	~EulerK();

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

//----------------------BilinearForm----------------------------------------------
 class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int entry_i, int entry_j, double kappa) : MatrixFormVol<double>(entry_i,entry_j),entry_i(entry_i), entry_j(entry_j), kappa(kappa) {}


    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;
		int entry_j; 
	int entry_i;
	double kappa;
	};

};





class EulerS : public WeakForm<double>
{
public:

  EulerS(double kappa,MeshFunctionSharedPtr<double>  rho_ext, MeshFunctionSharedPtr<double>  v1_ext, MeshFunctionSharedPtr<double>  v2_ext, MeshFunctionSharedPtr<double>  energy_ext, 
MeshFunctionSharedPtr<double>  prev_density, MeshFunctionSharedPtr<double>  prev_density_vel_x,  MeshFunctionSharedPtr<double>  prev_density_vel_y, MeshFunctionSharedPtr<double>  prev_energy, 
bool mirror_condition= true , int num_of_equations = 4);

	~EulerS();

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


//---Bilinear form bdry
 class EulerBoundaryBilinearForm: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm(double kappa,int entry_i, int entry_j) : MatrixFormSurf<double>(entry_i,entry_j),entry_i(entry_i), entry_j(entry_j), kappa(kappa) {}
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

	int entry_j; 
	int entry_i;
	double kappa;
	};



//--------------linearform boundary
 class EulerBoundaryLinearform: public VectorFormSurf<double>
  {
  public:
    EulerBoundaryLinearform(double kappa, int entry_i) : VectorFormSurf<double>(entry_i), kappa(kappa), entry_i(entry_i)  {}

     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double kappa;
		int entry_i;
	};



};








#endif

