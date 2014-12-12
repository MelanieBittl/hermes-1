#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "euler_util.h"
#include "euler_flux.h"
#include "hermes2d.h"
 #define PI (3.141592653589793)   



using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;

//---------------Massmatrix Part---------------
class EulerEquationsWeakForm_Mass : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_Mass(int num_of_equations = 8);

    WeakForm<double>* clone() const;

int num_of_equations;

};

//----------------Konvektion Part------------------------
class EulerK : public WeakForm<double>
{
public:

  EulerK(double gamma,
MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p, int num_of_equations = 8);

	~EulerK();

    WeakForm<double>* clone() const;


  // Members.
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	double gamma;
	

  MeshFunctionSharedPtr<double> prev_density_g;
  MeshFunctionSharedPtr<double> prev_density_vel_x_g;
  MeshFunctionSharedPtr<double> prev_density_vel_y_g;
  MeshFunctionSharedPtr<double> prev_energy_g;

  MeshFunctionSharedPtr<double> prev_density_p;
  MeshFunctionSharedPtr<double> prev_density_vel_x_p;
  MeshFunctionSharedPtr<double> prev_density_vel_y_p;
  MeshFunctionSharedPtr<double> prev_energy_p;


protected:

//----------------------BilinearForm Konvektion----------------------------------------------
 class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int entry_i, int entry_j, double gamma, bool particle = false) : MatrixFormVol<double>(entry_i,entry_j),entry_i(entry_i), entry_j(entry_j), gamma(gamma), particle(particle) {}


    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;
		int entry_j; 
	int entry_i;
	double gamma;
bool particle;
	};



 class  EulerEquationsLinearForm: public VectorFormVol<double>
  {
  public:
   EulerEquationsLinearForm(int entry_i, double gamma,bool particle=false) : VectorFormVol<double>(entry_i), entry_i(entry_i), gamma(gamma),particle(particle) {}

     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormVol<double>* clone() const;

	int entry_i;
	double gamma;
bool particle;

	};

};


//----------------SourceTerm Part------------------------
class EulerSource : public WeakForm<double>
{
public:

  EulerSource(double particle_density, double d, double c_vg,double c_vp, double c_pg, double Pr, double mu, MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p, int num_of_equations = 8);

	~EulerSource();

    WeakForm<double>* clone() const;

	double particle_density;

		double d;		
		double c_vg;
		double c_vp;
		double c_pg;
		double Pr;		
		double mu;


  MeshFunctionSharedPtr<double> prev_density_g;
  MeshFunctionSharedPtr<double> prev_density_vel_x_g;
  MeshFunctionSharedPtr<double> prev_density_vel_y_g;
  MeshFunctionSharedPtr<double> prev_energy_g;

  MeshFunctionSharedPtr<double> prev_density_p;
  MeshFunctionSharedPtr<double> prev_density_vel_x_p;
  MeshFunctionSharedPtr<double> prev_density_vel_y_p;
  MeshFunctionSharedPtr<double> prev_energy_p;


protected:

//----------------------BilinearForm----------------------------------------------
 class EulerSourceBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerSourceBilinearForm(int entry_i, int entry_j) : MatrixFormVol<double>(entry_i,entry_j),entry_i(entry_i), entry_j(entry_j) {}


    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;
		int entry_j; 
	int entry_i;

	};



 class  EulerSourceLinearForm: public VectorFormVol<double>
  {
  public:
   EulerSourceLinearForm(int entry_i) : VectorFormVol<double>(entry_i), entry_i(entry_i) {}

     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormVol<double>* clone() const;

	int entry_i;

	};

};



//------------Boundary+Penalty -----------------

class EulerBoundary : public WeakForm<double>
{
public:

  EulerBoundary(double gamma,double sigma,double particle_density,MeshFunctionSharedPtr<double>  rho_ext_g, MeshFunctionSharedPtr<double>  v1_ext_g, MeshFunctionSharedPtr<double>  v2_ext_g, MeshFunctionSharedPtr<double>  energy_ext_g,MeshFunctionSharedPtr<double>  rho_ext_p, MeshFunctionSharedPtr<double>  v1_ext_p, MeshFunctionSharedPtr<double>  v2_ext_p, MeshFunctionSharedPtr<double>  energy_ext_p, 
MeshFunctionSharedPtr<double>  prev_density_g, MeshFunctionSharedPtr<double>  prev_density_vel_x_g,  MeshFunctionSharedPtr<double>  prev_density_vel_y_g, MeshFunctionSharedPtr<double>  prev_energy_g,MeshFunctionSharedPtr<double>  prev_density_p, MeshFunctionSharedPtr<double>  prev_density_vel_x_p,  MeshFunctionSharedPtr<double>  prev_density_vel_y_p, MeshFunctionSharedPtr<double>  prev_energy_p,double eps = 1e-8,  int num_of_equations = 8);

	~EulerBoundary();

    WeakForm<double>* clone() const;



	   // Members.
	double sigma;//Penalty parameter
	double eps;
double particle_density;
  EulerFluxes* euler_fluxes;
	RiemannInvariants* riemann_invariants;
	double gamma;

  MeshFunctionSharedPtr<double> prev_density_g;
  MeshFunctionSharedPtr<double> prev_density_vel_x_g;
  MeshFunctionSharedPtr<double> prev_density_vel_y_g;
  MeshFunctionSharedPtr<double> prev_energy_g;

  MeshFunctionSharedPtr<double> prev_density_p;
  MeshFunctionSharedPtr<double> prev_density_vel_x_p;
  MeshFunctionSharedPtr<double> prev_density_vel_y_p;
  MeshFunctionSharedPtr<double> prev_energy_p;

  MeshFunctionSharedPtr<double> rho_ext_g;
  MeshFunctionSharedPtr<double> v1_ext_g;
  MeshFunctionSharedPtr<double> v2_ext_g;
  MeshFunctionSharedPtr<double> energy_ext_g;

  MeshFunctionSharedPtr<double> rho_ext_p;
  MeshFunctionSharedPtr<double> v1_ext_p;
  MeshFunctionSharedPtr<double> v2_ext_p;
  MeshFunctionSharedPtr<double> energy_ext_p;

protected:


//---Bilinear form bdry
 class EulerBoundaryBilinearForm: public MatrixFormSurf<double>
  {
  public:
    EulerBoundaryBilinearForm(double gamma,int entry_i, int entry_j, bool particle = false) : MatrixFormSurf<double>(entry_i,entry_j),entry_i(entry_i), entry_j(entry_j), gamma(gamma), particle(particle) {}
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

	int entry_j; 
	int entry_i;
	double gamma;
	bool particle;
	};



//--------------linearform boundary
 class EulerBoundaryLinearform: public VectorFormSurf<double>
  {
  public:
    EulerBoundaryLinearform(double gamma, int entry_i, bool particle = false) : VectorFormSurf<double>(entry_i), gamma(gamma), entry_i(entry_i), particle(particle)  {}

     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

		double gamma;
		int entry_i;
		bool particle;
	};
	
	
	//---Bilinearform penalty
 class PenaltyBilinearForm: public MatrixFormSurf<double>
  {
  public:
    PenaltyBilinearForm(double sigma,double eps,int entry_i, int entry_j, bool particle = false) : MatrixFormSurf<double>(entry_i,entry_j),eps(eps),entry_i(entry_i), entry_j(entry_j), sigma(sigma), particle(particle) {}
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormSurf<double>* clone() const;

	int entry_j; 
	int entry_i;
	double sigma;
	double eps;
	double particle;
	};

//--------------linearform penalty
 class PenaltyLinearForm: public VectorFormSurf<double>
  {
  public:
    PenaltyLinearForm(double sigma, int entry_i, bool particle = false) : VectorFormSurf<double>(entry_i), sigma(sigma) , entry_i(entry_i), particle(particle)  {}

     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    VectorFormSurf<double>* clone() const;

	double sigma;
		int entry_i;
	double particle;
	};



};










#endif

