#ifndef __H2D_MIXINS_H
#define __H2D_MIXINS_H
#include "global.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class Space;
    namespace Mixins
    {
      template<typename Scalar>
      class HERMES_API SettableSpaces
      {
      public:
        /// Sets new spaces for the instance.
        virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces) = 0;
        virtual void set_space(const Space<Scalar>* space) = 0;
        /// Get all spaces as a Hermes::vector.
        virtual Hermes::vector<const Space<Scalar>*> get_spaces() const = 0;
        virtual const Space<Scalar>* get_space(int n) const;
      };

      class HERMES_API StateQueryable
      {
      public:
        /// Ask if the instance is fine.
        virtual bool isOkay() const = 0;
      };

      template<typename Scalar>
      class HERMES_API MatrixRhsOutput
      {
      public:
        /// Constructor.
        /// Sets defaults (see individual set methods for values of those).
        MatrixRhsOutput();
        /// Sets this instance to output the matrix in several first iterations.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the matrix will be saved.
        void output_matrix(int firstIterations = -1);
        /// Sets filename for the matrix
        /// Default: Matrix_'iteration number' with the ".m" extension in the case of matlab format.
        /// \param[in] name sets the main part of the name, i.e. replacement for "Matrix_" in the default name.
        void set_matrix_filename(std::string name);
        /// Sets varname for the matrix
        /// Default: "A".
        void set_matrix_varname(std::string name);
        /// Sets varname for the matrix
        /// Default: "DF_MATLAB_SPARSE - matlab file".
        void set_matrix_E_matrix_dump_format(EMatrixDumpFormat format);
        
        /// Sets this instance to output the rhs in several first iterations.
        /// \param[in] firstIterations Only during so many first iterations. Default: -1 meaning, that during all iterations, the rhs will be saved.
        void output_rhs(int firstIterations = -1);
        /// Sets filename for the rhs
        /// Default: Rhs_'iteration number' with the ".m" extension in the case of matlab format.
        /// \param[in] name sets the main part of the name, i.e. replacement for "Rhs_" in the default name.
        void set_rhs_filename(std::string name);
        /// Sets varname for the rhs
        /// Default: "b".
        void set_rhs_varname(std::string name);
        /// Sets varname for the rhs
        /// Default: "DF_MATLAB_SPARSE - matlab file".
        void set_rhs_E_matrix_dump_format(EMatrixDumpFormat format);
        
      protected:
        bool output_matrixOn;
        int output_matrixIterations;
        std::string matrixFilename;
        std::string matrixVarname;
        EMatrixDumpFormat matrixFormat;

        bool output_rhsOn;
        int output_rhsIterations;
        std::string RhsFilename;
        std::string RhsVarname;
        EMatrixDumpFormat RhsFormat;
      };
    }
  }
}
#endif