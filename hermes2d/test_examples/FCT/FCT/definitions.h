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
      : MatrixFormVol<double>(i, j), time_step(time_step) { };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;

    // Members.  
    double time_step;
};


class VectorFormVolMass : public VectorFormVol<double>
{
public:
	VectorFormVolMass(int i, double  time_step) : VectorFormVol<double>(i), time_step(time_step) { };

	template<typename Real, typename Scalar>
	Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

	virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

	virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

   VectorFormVol<double>* clone() const;

	// Members.  
	double time_step;

};

class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix(double time_step,MeshFunctionSharedPtr<double> sln_prev_time);
//	~CustomWeakFormMassmatrix();
};

//---------------Konvektion-----------
class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolConvection(int i, int j) 
    : MatrixFormVol<double>(i, j) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, Func<Ord> **ext) const;  

    MatrixFormVol<double>* clone() const;

};


class VectorFormVolConvection : public VectorFormVol<double>
{
public:
  VectorFormVolConvection(int i) : VectorFormVol<double>(i){ }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

   VectorFormVol<double>* clone() const;
};


class CustomWeakFormConvection : public WeakForm<double>    //Konvektion
{
public:
  CustomWeakFormConvection(MeshFunctionSharedPtr<double> sln_prev_time);
	//~CustomWeakFormConvection();  
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
  GradientReconstruction_1(MeshFunctionSharedPtr<double> sln);
	//~GradientReconstruction_1();
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
  GradientReconstruction_2(MeshFunctionSharedPtr<double> sln);
	//~GradientReconstruction_2();
};
//----------------Filter-------------
    class AbsDifffilter : public DiffFilter<double>
    {
    public:
      AbsDifffilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions,  Hermes::vector<int> items = *(new Hermes::vector<int>)):DiffFilter(solutions,items){};
      virtual MeshFunction<double>* clone() const
		{
		  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
		  Hermes::vector<int> items;
		  for(int i = 0; i < this->num; i++)
		  {
		    slns.push_back(this->sln[i]->clone());
		    items.push_back(this->item[i]);
		  }
		  AbsDifffilter* filter = new AbsDifffilter(slns, items);		
		  return filter;
		}

virtual ~AbsDifffilter(){};
    protected:
      virtual void filter_fn(int n, Hermes::vector<double*> values, double* result)
    	{
      for (int i = 0; i < n; i++) result[i] = std::abs(values.at(0)[i] - values.at(1)[i]);
    	}
    };

//------------------- Initial condition ----------------

class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
   ~CustomInitialCondition(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;

   MeshFunction<double>* clone() const;
};



#endif

