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
                             double* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,bool* fct = nullptr, CSCMatrix<double>*  mat  =  nullptr);


	void project_lumped(std::vector<SpaceSharedPtr<double> >  spaces, std::vector<MeshFunctionSharedPtr<double> > source_meshfns,
          double* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,bool* fct = nullptr,CSCMatrix<double>*  mat  =  nullptr);


protected:
 	void project_internal( SpaceSharedPtr<double> space, WeakForm<double>* wf, double* target_vec,
                               MatrixSolverType matrix_solver,bool* fct = nullptr, CSCMatrix<double>*  mat =  nullptr);


};


#endif
