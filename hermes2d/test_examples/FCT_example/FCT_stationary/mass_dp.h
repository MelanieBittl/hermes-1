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

class HERMES_API Lumped_Flux : public DiscreteProblem<double>  
 { 
   public:
		      /// Constructor for one equation.
      Lumped_Flux(const WeakForm<double>* wf, const Space<double>* space, double* u_H, double* u_L) : DiscreteProblem<double>(wf, space), u_H(u_H),u_L(u_L){
			ref_ndof = space->get_num_dofs();
			P_plus = new double[ref_ndof];  P_minus = new double[ref_ndof];
			Q_plus = new double[ref_ndof]; Q_minus = new double[ref_ndof];	
			R_plus = new double[ref_ndof]; R_minus = new double[ref_ndof];
			 memset(P_plus, 0, ref_ndof*sizeof(double));memset(P_minus, 0, ref_ndof*sizeof(double));
			memset(Q_plus, 0, ref_ndof*sizeof(double));memset(Q_minus, 0, ref_ndof*sizeof(double));
			memset(R_plus, 0, ref_ndof*sizeof(double));memset(R_minus, 0, ref_ndof*sizeof(double));
			elem_flux = new double[4];	memset(elem_flux, 0, 4*sizeof(double));
		};

      /// Destuctor.
       ~Lumped_Flux(){};

      /// Assembling.
      /// General assembling procedure for nonlinear problems. coeff_vec is the
      /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems. The parameter add_dir_lift decides
      /// whether Dirichlet lift will be added while coeff_vec is converted into
      /// Solutions.
			void assemble(double* coeff_vec, SparseMatrix<double>* mat, Vector<double>* rhs = NULL,
        bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL);

      /// Assembling.
      /// Without the matrix.
      void assemble(double* coeff_vec, Vector<double>* rhs = NULL,
        bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      void assemble(SparseMatrix<double>* mat, Vector<double>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

		 double* elem_flux;		
		double* P_plus; double* P_minus; 
		double* Q_plus ; double * Q_minus ; 
		double* R_plus; double* R_minus; 
		int ref_ndof;
protected:
      /// Assemble volume matrix forms.
      void assemble_volume_matrix_forms(Stage<double>& stage,
        SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<double>*>& u_ext,
        int marker, Hermes::vector<AsmList<double>*>& al);
/// Assemble vector matrix forms.
       void assemble_volume_vector_forms(Stage<double>& stage,
      SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks,
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss,
      Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution<double>*>& u_ext,
      int marker, Hermes::vector<AsmList<double>*>& al);

    // Members.  
    double* u_H;
		double* u_L;


};

class HERMES_API Artificial_Diffusion : public DiscreteProblem<double>  
 { 
   public:
		      /// Constructor for one equation.
      Artificial_Diffusion (const WeakForm<double>* wf, const Space<double>* space) : DiscreteProblem<double>(wf, space){
		};

      /// Destuctor.
       ~Artificial_Diffusion(){};

protected:
      /// Assemble volume matrix forms.
      void assemble_volume_matrix_forms(Stage<double>& stage,
        SparseMatrix<double>* mat, Vector<double>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<double>*>& u_ext,
        int marker, Hermes::vector<AsmList<double>*>& al);


};






#endif

