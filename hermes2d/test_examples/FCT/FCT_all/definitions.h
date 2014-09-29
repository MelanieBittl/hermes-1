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

class CustomMatrixFormVolDiffusion : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolDiffusion(int i, int j) 
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

class CustomMatrixFormVolReaction : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolReaction(int i, int j) 
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


class CustomWeakFormConvection : public WeakForm<double>    //Konvektion
{
public:
  CustomWeakFormConvection(bool conv, bool diff,bool reac);
};


//--------surface----------------

class CustomWeakFormSurface : public WeakForm<double>    //Surface
{
public:
  CustomWeakFormSurface(MeshFunctionSharedPtr<double> sln_prev_time);

	
private:

  class CustomMatrixFormSurface : public MatrixFormSurf<double>
  {
  public:
    CustomMatrixFormSurface(int i, int j) : MatrixFormSurf<double>(i, j) {};


    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormSurf<double>* clone() const;

  };

  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface(int i) : VectorFormSurf<double>(i) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;

  };
  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;
};
//---------RHS force---------
class CustomWeakFormRHS : public WeakForm<double>
{
public:
  CustomWeakFormRHS(MeshFunctionSharedPtr<double> sln_exact,MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,double time_step,double theta);
  WeakForm<double>* clone() const;
	
private:

  class RHS : public VectorFormVol<double>
  {
  public:
    RHS(int i, double time_step, double theta) : VectorFormVol<double>(i), time_step(time_step), theta(theta)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;
double time_step;
double theta;

  };

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
  CustomInitialCondition(MeshSharedPtr mesh, double time) : ExactSolutionScalar<double>(mesh) , time(time)
{
		
};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;

   MeshFunction<double>* clone() const;
	void set_time(double t){time = t;}

protected:
double time;
};



#endif

