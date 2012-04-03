#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "hermes2d.h"
 #define PI (3.141592653589793)   



using namespace Hermes;
using namespace Hermes::Hermes2D;


//---------------Massematrix-----------

class CustomMatrixFormVolMassmatrix : public MatrixFormVol<double>   

{
  public:
    // This weak form is custom since it contains a nonlinearity in the diffusion term.
    CustomMatrixFormVolMassmatrix(int i, int j, double time_step) 
      : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step) { };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();

    // Members.  
    double time_step;
};


class VectorFormVolMass : public VectorFormVol<double>
{
public:
	VectorFormVolMass(int i, double  time_step) : VectorFormVol<double>(i), time_step(time_step) { };

	template<typename Real, typename Scalar>
	Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

	virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();

	// Members.  
	double time_step;

};

class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix(double time_step,Solution<double>* sln_prev_time);
	~CustomWeakFormMassmatrix();
};

//---------------Konvektion-----------
class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolConvection(int i, int j) 
    : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;  

    MatrixFormVol<double>* clone();

};


class VectorFormVolConvection : public VectorFormVol<double>
{
public:
  VectorFormVolConvection(int i) : VectorFormVol<double>(i){ }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();
};


class CustomWeakFormConvection : public WeakForm<double>    //Konvektion
{
public:
  CustomWeakFormConvection(Solution<double>* sln_prev_time);
	~CustomWeakFormConvection();  
};


//------Matrix & Vektorform for higher Order solution----------------

 class ConvectionMatForm : public MatrixFormVol<double>
  {
  public:
    // This weak form is custom since it contains a nonlinearity in the diffusion term.
    ConvectionMatForm(int i, int j, double time_step) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const ;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const ;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();

    // Members.
  
    double time_step;
  };

  // This form (residual) is custom since it contains a nonlinear term.
 class VectorConvection : public VectorFormVol<double>
  {
  public:
    VectorConvection(int i,double time_step) : VectorFormVol<double>(i), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double>  *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();

    // Members.
    
    double time_step;
  };

class ConvectionForm : public WeakForm<double> 
{
public:
  ConvectionForm( double time_step, Solution<double>* sln_prev_time) ;
	~ConvectionForm();
  };






//-----------------Residual for error estimator--------------------------------
class ResidualMatForm : public MatrixFormVol<double>
{
public:
  ResidualMatForm(int i, int j, double time_step) 
    : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
  

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;    

    MatrixFormVol<double>* clone();

  // Members.  
  double time_step;
};


class VectorResidual : public VectorFormVol<double>
{
public:
  VectorResidual(int i,double time_step) : VectorFormVol<double>(i), time_step(time_step) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const ;
 

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();

  // Members.    
  double time_step;
};
class ResidualForm : public WeakForm<double>
{
public:
  ResidualForm( double time_step, Solution<double>* sln_prev_time, Solution<double>* ref);
	~ResidualForm();
};

//---------Gradient Reconstruction---------------

class GradientReconstructionMatForm_1 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_1(int i, int j) 
    : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
  

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;  
  
    MatrixFormVol<double>* clone();

};


class GradientReconstructionVectorForm_1 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_1(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const ;
 

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();

};
class GradientReconstruction_1 : public WeakForm<double>
{
public:
  GradientReconstruction_1(Solution<double>* sln);
	~GradientReconstruction_1();
};


class GradientReconstructionMatForm_2 : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm_2(int i, int j) 
    : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM){ }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;
  

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;    

    MatrixFormVol<double>* clone();
};


class GradientReconstructionVectorForm_2 : public VectorFormVol<double>
{
public:
  GradientReconstructionVectorForm_2(int i) : VectorFormVol<double>(i) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const ;
 

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

   VectorFormVol<double>* clone();

};
class GradientReconstruction_2 : public WeakForm<double>
{
public:
  GradientReconstruction_2(Solution<double>* sln);
	~GradientReconstruction_2();
};

//------------------- Initial condition ----------------

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(Ord x, Ord y) const ;

   MeshFunction<double>* clone();
};



#endif

