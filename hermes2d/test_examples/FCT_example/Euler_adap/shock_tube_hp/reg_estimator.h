#ifndef __REGESTIMATOR_H
#define __REGESTIMATOR_H


#include "hermes2d.h"
#include <list>


using namespace Hermes;
using namespace Hermes::Hermes2D;

//---------Gradient Reconstruction---------------
class GradientReconstructionMatForm : public MatrixFormVol<double>
{
public:
  GradientReconstructionMatForm(int i, int j) 
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


//----------------Regularity_Estimator Class


class Regularity_Estimator
{
public:

	Regularity_Estimator(double epsilon);	
	~Regularity_Estimator();
	void free();
	
	void set_space(Space<double>* new_space);
	
	int* get_smooth_elems(Space<double>* new_space,double* coeff_vec,UMFPackMatrix<double> * mass_matrix=NULL);
	int* get_smooth_dofs(Space<double>* new_space,double* coeff_vec,UMFPackMatrix<double> * mass_matrix=NULL);

	int* get_smooth_elems(Space<double>* new_space,Solution<double>* new_sln,UMFPackMatrix<double> * mass_matrix=NULL);
	int* get_smooth_dofs(Space<double>* new_space,Solution<double>* new_sln,UMFPackMatrix<double> * mass_matrix=NULL);


protected:
	double linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln);
	double linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln);
	double linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c,Solution<double>* sln);

	void smoothness_indicator(UMFPackMatrix<double> * mass_matrix=NULL);

	

	double epsilon;
	Solution<double>* sln;
	Solution<double>* R_h_1;
	Solution<double>* R_h_2;
	GradientReconstruction_1* grad_1;
	GradientReconstruction_2* grad_2;
	Space<double>* space;
	UMFPackVector<double> * rhs_1;
	UMFPackVector<double> * rhs_2;

	int* smooth_elem_patch;
	int* smooth_dof;
		AsmList<double>* al;

};

#endif
