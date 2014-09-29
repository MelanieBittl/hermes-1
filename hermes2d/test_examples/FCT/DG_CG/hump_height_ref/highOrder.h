#ifndef __HIGHORDER_H
#define __HIGHORDER_H

#include "hermes2d.h"
using namespace Hermes::Solvers;
class High_Order{
public:
	//Constructor
	High_Order(double theta);

	//Destructor
	~High_Order();

	void assemble_High_Order(CSCMatrix<double> * conv_matrix, CSCMatrix<double> * mass_matrix,CSCMatrix<double> * dg_surface_matrix,CSCMatrix<double> * reac_matrix);
	double* solve_High_Order(double* u_n,SimpleVector<double> * rhs_f);


protected:

	CSCMatrix<double> * high_matrix;  
	CSCMatrix<double> * highmat_rhs;
	double* u_H;
	double theta;

};






#endif

