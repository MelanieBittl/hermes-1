#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "hermes2d.h"
 #define PI (3.141592653589793) 

using namespace Hermes;
using namespace Hermes::Hermes2D;






//---------------mass-matrix/tau-----------

class CustomMatrixFormVolMassmatrix : public MatrixFormVol<double>   

{
  public:
    CustomMatrixFormVolMassmatrix(int i, int j) 
      : MatrixFormVol<double>(i, j) { };

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const ;

};





class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix();
	
};




//--------------------Error_calculation--------------

class CustomNormFormVol : public NormFormVol<double>
{
public:
  CustomNormFormVol(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class StreamlineDiffusionNorm : public NormFormVol<double>
{
public:
  StreamlineDiffusionNorm(int i, int j,  MeshSharedPtr mesh);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;

protected:
MeshSharedPtr mesh;

};

class CustomNormFormSurf : public NormFormSurf<double>
{
public:
  CustomNormFormSurf(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class CustomNormFormDG : public NormFormDG<double>
{
public:
  CustomNormFormDG(int i, int j);

  virtual double value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const;
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
  CustomInitialCondition(MeshSharedPtr mesh, bool all) : ExactSolutionScalar<double>(mesh), all(all) {};


  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

 virtual Ord ord(double x, double y)  const ;

   MeshFunction<double>* clone() const;
	bool all;
};




#endif

