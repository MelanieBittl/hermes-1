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




//----------------Bilinearform for K

 class EulerEquationsBilinearFormDensity : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensity(int j, EulerFluxes* euler_fluxes) : MatrixFormVol<double>(0,j),euler_fluxes(euler_fluxes), j(j) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const ;    
     

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const; 


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const; 

    MatrixFormVol<double>* clone() const;
		int j;
EulerFluxes* euler_fluxes;
	};

 class EulerEquationsBilinearFormDensityVelX : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensityVelX(double kappa, int j, EulerFluxes* euler_fluxes)
      : MatrixFormVol<double>(1,j), kappa(kappa),euler_fluxes(euler_fluxes), j(j) {}

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
EulerFluxes* euler_fluxes;
  };

  class EulerEquationsBilinearFormDensityVelY : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormDensityVelY(double kappa,int j, EulerFluxes* euler_fluxes) 
      : MatrixFormVol<double>(2,j), kappa(kappa),euler_fluxes(euler_fluxes), j(j) {}

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
EulerFluxes* euler_fluxes;
  };


 class EulerEquationsBilinearFormEnergy : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormEnergy(double kappa, int j, EulerFluxes* euler_fluxes) 
      : MatrixFormVol<double>(3,j), kappa(kappa),	euler_fluxes(euler_fluxes), j(j) {}

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
EulerFluxes* euler_fluxes;
  };





//---------------WeakForms---------------
class EulerEquationsWeakForm_Mass : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_Mass(int num_of_equations = 4);

    WeakForm<double>* clone() const;

 int num_of_equations;
};

class EulerEquationsWeakForm_K : public WeakForm<double>
{
public:
  EulerEquationsWeakForm_K(double kappa, Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns, int num_of_equations = 4);
    WeakForm<double>* clone() const;
 ~EulerEquationsWeakForm_K();

  // Members.
  EulerFluxes* euler_fluxes;
	double kappa;
Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns;

};




//---------Gradient Reconstruction---------------

class GradientReconstructionMatForm_1 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_1(int i, int j)
    : MatrixFormVol<double>(i, j){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;
  

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
               Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, Func<Ord> **ext) const;
  
    MatrixFormVol<double>* clone() const;

};


class GradientReconstructionVectorForm_1 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_1(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const ;
 

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

   VectorFormVol<double>* clone() const;

};
class GradientReconstruction_1 : public WeakForm<double>
{
public:
  GradientReconstruction_1( MeshFunctionSharedPtr<double> sln);
      WeakForm<double>* clone() const;
~GradientReconstruction_1();
};


class GradientReconstructionMatForm_2 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_2(int i, int j)
    : MatrixFormVol<double>(i, j){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;
  

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
               Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;
};


class GradientReconstructionVectorForm_2 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_2(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const ;
 

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

   VectorFormVol<double>* clone() const;

};
class GradientReconstruction_2 : public WeakForm<double>
{
public:
  GradientReconstruction_2( MeshFunctionSharedPtr<double> sln);
      WeakForm<double>* clone() const;
~GradientReconstruction_2();
};



#endif

