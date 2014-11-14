#ifndef __LOWORDER_H
#define __LOWORDER_H

#include "hermes2d.h"

using namespace Hermes::Solvers;
class Low_Order{

public:
	//Constructor
	Low_Order(double theta);

	//Destructor
	~Low_Order();

	void assemble_Low_Order(CSCMatrix<double> * conv_matrix,CSCMatrix<double>* diffusion, CSCMatrix<double> * lumped_matrix,CSCMatrix<double> * surface_matrix,CSCMatrix<double> * dg_matrix);
	double* solve_Low_Order(CSCMatrix<double> * lumped_matrix, double* u_n,SimpleVector<double> * surf, double time_step);
	double* explicit_Correction(double* flux_correction);


protected:

	CSCMatrix<double> * low_matrix;  
	CSCMatrix<double> * lowmat_rhs;
	double* u_L;
	double* rhs;
	double* u_new;	
	double theta;


};
























#endif

