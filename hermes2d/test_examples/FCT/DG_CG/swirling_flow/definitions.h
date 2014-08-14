#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "hermes2d.h"
 #define PI (3.141592653589793) 

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;


class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(MeshFunctionSharedPtr<double> sln_prev_time,MeshSharedPtr mesh,double final_time,double time_step,double theta, bool all = false, bool DG = true, bool SD =false, bool right_hand_side = true);
  WeakForm<double>* clone() const;
double final_time;
	
private:

class CustomMatrixFormVolConvection : public MatrixFormVol<double>   
{
public:
  
  CustomMatrixFormVolConvection(int i, int j, double time_step, double theta) 
    : MatrixFormVol<double>(i, j), time_step(time_step), theta(theta) { }

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
               Func<double> *v, Geom<double> *e, Func<double>  **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
          Geom<Ord> *e, Func<Ord>  **ext) const;  
    MatrixFormVol<double>* clone() const;
double time_step;
double theta;

};


  class Streamline : public MatrixFormVol<double>
  {
  public:
    Streamline(int i, int j,MeshSharedPtr mesh) : MatrixFormVol<double>(i, j), mesh(mesh) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;
		MeshSharedPtr mesh;
  };



  class CustomMatrixFormSurface : public MatrixFormSurf<double>
  {
  public:
    CustomMatrixFormSurface(int i, int j) : MatrixFormSurf<double>(i, j) {};


    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    MatrixFormSurf<double>* clone() const;

  };

  class CustomMatrixFormInterface : public MatrixFormDG<double>
  {
  public:
    CustomMatrixFormInterface(int i, int j, double theta ) : MatrixFormDG<double>(i, j) , theta(theta) 
    {
    };

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, DiscontinuousFunc<Ord> *u, DiscontinuousFunc<Ord> *v, Geom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    MatrixFormDG<double>* clone() const;
double theta;

  };

  class CustomVectorFormInterface : public VectorFormDG<double>
  {
  public:
			CustomVectorFormInterface(int i, double theta) : VectorFormDG<double>(i), theta(theta) 
    {
    };

      virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
        Geom<double> *e, DiscontinuousFunc<double> **ext) const;

      virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, Geom<Ord> *e,
        DiscontinuousFunc<Ord> **ext) const;

      virtual VectorFormDG<double>* clone() const;


double theta;

  };


  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface(int i) : VectorFormSurf<double>(i) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;

  };


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
  
  double calculate_a_dot_v(double x, double y, double vx, double vy) const;

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const;

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;
  
    double get_final_time(){ return final_time;}


};


//---------------Mass matrix

class  CustomWeakFormMassmatrix  : public WeakForm<double>     
{
public:
  CustomWeakFormMassmatrix();
    WeakForm<double>* clone() const;
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
  StreamlineDiffusionNorm(int i, int j,  MeshSharedPtr mesh, double final_time, double time);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;

protected:
MeshSharedPtr mesh;
double time;
double final_time;

};

class CustomNormFormSurf : public NormFormSurf<double>
{
public:
  CustomNormFormSurf(int i, int j,double final_time, double time);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
protected:
	double time;
double final_time;
};

class CustomNormFormDG : public NormFormDG<double>
{
public:
  CustomNormFormDG(int i, int j,double final_time, double time);

  virtual double value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const;
  protected:
	double time;
double final_time;
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

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const;

   MeshFunction<double>* clone() const;
};




#endif

