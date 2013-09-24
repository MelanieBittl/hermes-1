// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
/*! \file umfpack_solver.cpp
\brief UMFPACK solver interface.
*/
#include "config.h"
#ifdef WITH_UMFPACK
#include "umfpack_solver.h"
#include "common.h"

#define umfpack_real_symbolic umfpack_di_symbolic
#define umfpack_real_numeric umfpack_di_numeric
#define umfpack_real_solve umfpack_di_solve

#define umfpack_complex_symbolic umfpack_zi_symbolic
#define umfpack_complex_numeric umfpack_zi_numeric
#define umfpack_complex_solve umfpack_zi_solve

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    void UMFPackLinearMatrixSolver<Scalar>::set_output_level(double level)
    {
      Control[UMFPACK_PRL] = level;
    }

    template<typename Scalar>
    UMFPackLinearMatrixSolver<Scalar>::UMFPackLinearMatrixSolver(CSCMatrix<Scalar> *m, SimpleVector<Scalar> *rhs)
      : DirectSolver<Scalar>(HERMES_CREATE_STRUCTURE_FROM_SCRATCH), m(m), rhs(rhs), symbolic(NULL), numeric(NULL)
    {
      umfpack_di_defaults(Control);
    }

    template<typename Scalar>
    UMFPackLinearMatrixSolver<Scalar>::~UMFPackLinearMatrixSolver()
    {
      free_factorization_data();
    }

    template<typename Scalar>
    int UMFPackLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return m->get_size();
    }

    template<>
    bool UMFPackLinearMatrixSolver<double>::setup_factorization()
    {
      // Perform both factorization phases for the first time.
      if(reuse_scheme != HERMES_CREATE_STRUCTURE_FROM_SCRATCH && symbolic == NULL && numeric == NULL)
        reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH;
      else
        reuse_scheme = reuse_scheme;

      int status;
      switch(reuse_scheme)
      {
      case HERMES_CREATE_STRUCTURE_FROM_SCRATCH:
        if(symbolic != NULL)
        {
          umfpack_di_free_symbolic(&symbolic);
          memset(Info, 0, 90 * sizeof(double));
        }

        // Factorizing symbolically.
        status = umfpack_real_symbolic(m->get_size(), m->get_size(), m->get_Ap(), m->get_Ai(), m->get_Ax(), &symbolic, Control, Info);
        if(status != UMFPACK_OK)
        {
          if(symbolic)
            umfpack_di_free_symbolic(&symbolic);
          throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK symbolic factorization", status));
        }

      case HERMES_REUSE_MATRIX_REORDERING:
      case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
        if(numeric != NULL)
        {
          umfpack_di_free_numeric(&numeric);
          memset(Info + 0, 0, 90 * sizeof(double));
        }

        // Factorizing numerically.
        status = umfpack_real_numeric(m->get_Ap(), m->get_Ai(), m->get_Ax(), symbolic, &numeric, Control, Info);
        if(status != UMFPACK_OK)
        {
          if(numeric)
            umfpack_di_free_numeric(&numeric);
          throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK numeric factorization", status));
        }
        else
          umfpack_di_report_info (Control, Info);
      }

      return true;
    }

    template<>
    bool UMFPackLinearMatrixSolver<std::complex<double> >::setup_factorization()
    {
      // Perform both factorization phases for the first time.
      int eff_fact_scheme;
      if(reuse_scheme != HERMES_CREATE_STRUCTURE_FROM_SCRATCH && symbolic == NULL && numeric == NULL)
        eff_fact_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH;
      else
        eff_fact_scheme = reuse_scheme;

      int status;
      switch(eff_fact_scheme)
      {
      case HERMES_CREATE_STRUCTURE_FROM_SCRATCH:
        if(symbolic != NULL)
          umfpack_zi_free_symbolic(&symbolic);

        status = umfpack_complex_symbolic(m->get_size(), m->get_size(), m->get_Ap(), m->get_Ai(), (double *)m->get_Ax(), NULL, &symbolic, NULL, NULL);
        if(status != UMFPACK_OK)
        {
          if(symbolic)
            umfpack_zi_free_symbolic(&symbolic);
          throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK symbolic factorization", status));
        }

      case HERMES_REUSE_MATRIX_REORDERING:
      case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
        if(numeric != NULL)
          umfpack_zi_free_numeric(&numeric);

        status = umfpack_complex_numeric(m->get_Ap(), m->get_Ai(), (double *) m->get_Ax(), NULL, symbolic, &numeric, NULL, NULL);
        if(status != UMFPACK_OK)
        {
          if(numeric)
            umfpack_zi_free_numeric(&numeric);
          throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK numeric factorization", status));
        }
      }

      return true;
    }

    template<>
    void UMFPackLinearMatrixSolver<double>::free_factorization_data()
    {
      if(symbolic != NULL) umfpack_di_free_symbolic(&symbolic);
      symbolic = NULL;
      if(numeric != NULL) umfpack_di_free_numeric(&numeric);
      numeric = NULL;
    }

    template<>
    void UMFPackLinearMatrixSolver<std::complex<double> >::free_factorization_data()
    {
      if(symbolic != NULL) umfpack_zi_free_symbolic(&symbolic);
      symbolic = NULL;
      if(numeric != NULL) umfpack_zi_free_numeric(&numeric);
      numeric = NULL;
    }

    template<>
    void UMFPackLinearMatrixSolver<double>::solve()
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->get_size() == rhs->get_size());

      this->tick();

      if( !setup_factorization() )
        throw Exceptions::LinearMatrixSolverException("LU factorization could not be completed.");

      if(sln != NULL)
        delete [] sln;

      sln = new double[m->get_size()];
      memset(sln, 0, m->get_size() * sizeof(double));
      int status = umfpack_real_solve(UMFPACK_A, m->get_Ap(), m->get_Ai(), m->get_Ax(), sln, rhs->v, numeric, NULL, NULL);
      if(status != UMFPACK_OK)
      {
        this->free_factorization_data();
        throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK solution", status));
      }

      this->tick();
    }

    template<>
    void UMFPackLinearMatrixSolver<std::complex<double> >::solve()
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->get_size() == rhs->get_size());

      this->tick();
      if( !setup_factorization() )
        this->warn("LU factorization could not be completed.");

      if(sln)
        delete [] sln;
      sln = new std::complex<double>[m->get_size()];

      memset(sln, 0, m->get_size() * sizeof(std::complex<double>));
      int status = umfpack_complex_solve(UMFPACK_A, m->get_Ap(), m->get_Ai(), (double *)m->get_Ax(), NULL, (double*) sln, NULL, (double *)rhs->v, NULL, numeric, NULL, NULL);
      if(status != UMFPACK_OK)
      {
        this->free_factorization_data();
        throw Exceptions::LinearMatrixSolverException(check_status("UMFPACK solution", status));
      }

      this->tick();
      time = this->accumulated();
    }

    template<typename Scalar>
    char* UMFPackLinearMatrixSolver<Scalar>::check_status(const char *fn_name, int status)
    {
      char* to_return = new char[100];

      switch (status)
      {
      case UMFPACK_OK: break;
      case UMFPACK_WARNING_singular_matrix: sprintf(to_return, "%s: UMFPACK_WARNING_singular_matrix!", fn_name); break;
      case UMFPACK_ERROR_out_of_memory: sprintf(to_return, "%s: UMFPACK_ERROR_out_of_memory!", fn_name); break;
      case UMFPACK_ERROR_argument_missing: sprintf(to_return, "%s: UMFPACK_ERROR_argument_missing", fn_name); break;
      case UMFPACK_ERROR_invalid_Symbolic_object: sprintf(to_return, "%s: UMFPACK_ERROR_invalid_Symbolic_object", fn_name); break;
      case UMFPACK_ERROR_invalid_Numeric_object: sprintf(to_return, "%s: UMFPACK_ERROR_invalid_Numeric_object", fn_name); break;
      case UMFPACK_ERROR_different_pattern: sprintf(to_return, "%s: UMFPACK_ERROR_different_pattern", fn_name); break;
      case UMFPACK_ERROR_invalid_system: sprintf(to_return, "%s: UMFPACK_ERROR_invalid_system", fn_name); break;
      case UMFPACK_ERROR_n_nonpositive: sprintf(to_return, "%s: UMFPACK_ERROR_n_nonpositive", fn_name); break;
      case UMFPACK_ERROR_invalid_matrix: sprintf(to_return, "%s: UMFPACK_ERROR_invalid_matrix", fn_name); break;
      case UMFPACK_ERROR_internal_error: sprintf(to_return, "%s: UMFPACK_ERROR_internal_error", fn_name); break;
      default: sprintf(to_return, "%s: unknown error (%d)", fn_name, status); break;
      }

      return to_return;
    }

    template class HERMES_API UMFPackLinearMatrixSolver<double>;
    template class HERMES_API UMFPackLinearMatrixSolver<std::complex<double> >;
  }
}
#endif
