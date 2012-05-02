#ifndef __H2D_MASS_DP_H
#define __H2D_MASS_DP_H
#include "hermes2d.h"
#include "discrete_problem.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

    class Hermes::Hermes2D::PrecalcShapeset;

    class HERMES_API Mass_DP : public DiscreteProblem<double>  
 {    public:

		      /// Constructor for one equation.
      Mass_DP(const WeakForm<double>* wf, const Space<double>* space) : DiscreteProblem<double>(wf, space){
		};

      /// Destuctor.
       ~Mass_DP(){};

protected:
      /// Assemble volume matrix forms.
      void assemble_volume_matrix_forms(Stage<double>& stage,
        SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<double>*>& u_ext,
        int marker, Hermes::vector<AsmList<double>*>& al);


};










#endif

