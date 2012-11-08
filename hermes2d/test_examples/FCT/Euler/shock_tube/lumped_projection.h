#ifndef __LUMPED_PROJECTION_H
#define __LUMPED_PROJECTION_H


#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;




class Lumped_Projection
{
public:
 void project_lumped( const Space<double>* space, MeshFunction<double>* source_meshfn,
                             double* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,UMFPackMatrix<double>*  mat  = NULL);


void project_lumped(Hermes::vector<const Space<double>*> spaces, Hermes::vector<MeshFunction<double>*> source_meshfns,
          double* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,UMFPackMatrix<double>*  mat  = NULL);



protected:
 void project_internal( const Space<double>* space, WeakForm<double>* wf, double* target_vec,
                               MatrixSolverType matrix_solver, UMFPackMatrix<double>*  mat = NULL);


  class ProjectionLumpedMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    ProjectionLumpedMatrixFormVol(int i, int j) : MatrixFormVol<double>(i, j)
    {
    
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                 Geom<double> *e, Func<double>  **ext) const
    {      
      
        return lumped_projection_biform<double, double>(n, wt, u_ext, u, v, e, ext);
     
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Hermes::Ord>  **ext) const
    {
          
        return lumped_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
     
    }
    MatrixFormVol<double>* clone() const
    {
      return new ProjectionLumpedMatrixFormVol(*this);
    }


  private:
  

		template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain lumped_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) 
        {
        
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i]);
          return result;
        }

  };

  // Residual.
  class ProjectionLumpedVectorFormVol : public VectorFormVol<double>
  {
  public:
    ProjectionLumpedVectorFormVol(int i) : VectorFormVol<double>(i)
    {

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                 Geom<double> *e, Func<double>  **ext) const
    {
     
      
        return lumped_projection_residual<double, double>(n, wt, u_ext, v, e, ext);
     
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Hermes::Ord>  **ext) const
    {
      
        return lumped_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
    
    }

        VectorFormVol<double>* clone() const
        {
          return new ProjectionLumpedVectorFormVol(*this);
        }

  private:

   template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain lumped_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
 
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
							result += wt[i] * (ext[0]->val[i]) * v->val[i];

          return result;
        }
   
  };
};


#endif
