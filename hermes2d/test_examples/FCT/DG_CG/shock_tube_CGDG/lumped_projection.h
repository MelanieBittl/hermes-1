#ifndef __LUMPED_PROJECTION_H
#define __LUMPED_PROJECTION_H

#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Solvers;

class Lumped_Projection
{
public:
 	void project_lumped( SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double>  source_meshfn,
                             double* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,CSCMatrix<double>*  mat  = NULL);


	void project_lumped(Hermes::vector<SpaceSharedPtr<double> >  spaces, Hermes::vector<MeshFunctionSharedPtr<double> > source_meshfns,
          double* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,CSCMatrix<double>*  mat  = NULL);


protected:
 	void project_internal( SpaceSharedPtr<double> space, WeakForm<double>* wf, double* target_vec,
                               MatrixSolverType matrix_solver, CSCMatrix<double>*  mat = NULL);


};


#endif
