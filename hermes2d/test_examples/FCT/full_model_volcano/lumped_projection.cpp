#include "lumped_projection.h"
#include "discrete_problem.h"


////template<typename double>
void Lumped_Projection::project_internal( SpaceSharedPtr<double> space,WeakFormSharedPtr<double>  wf, double* target_vec,
                             CSCMatrix<double>*  mat)
{
   if (space == nullptr) printf("this->space == nullptr in Lumped_Projection::project_internal().");
  // Get dimension of the space.
  int ndof = space->get_num_dofs();

  if (mat != nullptr) if (mat->get_size() != ndof) printf("matrix=%i, ndof=%i", mat->get_size(), ndof);

  // Initialize DiscreteProblem.
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(wf, space);
  CSCMatrix<double>* matrix = new CSCMatrix < double > ;
  SimpleVector<double>* rhs = new SimpleVector<double>(ndof);
  double* coeff_vec = nullptr;
  if (mat == nullptr) 		//=> masslumping
  {
    //M_L
    CSCMatrix<double>* lumped_matrix = new CSCMatrix < double > ;
    dp->assemble(matrix, rhs);
    int size = matrix->get_size();
    double* diag = new double[size];
    int nnz = matrix->get_nnz();
    int* row = new int[size];
    int* col = new int[size + 1];
    for (int i = 0; i < size; i++)
    {
      diag[i] = 0;
      row[i] = i;
      col[i] = i;
    }
    col[size] = size;

    for (int i = 0; i < nnz; i++)
      diag[matrix->get_Ai()[i]] += matrix->get_Ax()[i];
    lumped_matrix->create(size, size, col, row, diag);
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(lumped_matrix, rhs);
    try
    {
      solver->solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    }
    coeff_vec = solver->get_sln_vector();

    if (target_vec != nullptr)
      for (int i = 0; i < ndof; i++)
        target_vec[i] = coeff_vec[i];
    delete solver;
    delete lumped_matrix;
    delete[] diag;
    delete[] row;
    delete[] col;
  }
  else
  {
    dp->assemble(rhs);
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(mat, rhs);

    try
    {
      solver->solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    }

    coeff_vec = solver->get_sln_vector();

    if (target_vec != nullptr)
      for (int i = 0; i < ndof; i++) target_vec[i] = coeff_vec[i];
    delete solver;
  }

  delete matrix;
  delete rhs;
  delete dp;
  
}


  void Lumped_Projection::project_lumped(std::vector<SpaceSharedPtr<double> >  spaces, 
        std::vector<MeshFunctionSharedPtr<double> > source_meshfns,
        double* target_vec,CSCMatrix<double>*  mat)
    {
     
      int n = spaces.size();

      // Sanity checks.
      if (n != source_meshfns.size()) throw Exceptions::LengthException(1, 2, n, source_meshfns.size());
      if (target_vec == NULL) throw Exceptions::NullException(3);
  
      int start_index = 0;
      for (int i = 0; i < n; i++) 
			{     
				spaces[i]->assign_dofs();
          project_lumped(spaces[i], source_meshfns[i], target_vec + start_index, mat);
        start_index += spaces[i]->get_num_dofs();
			}
		Space<double>::assign_dofs(spaces);
    }



void Lumped_Projection::project_lumped(const  SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double> source_meshfn,
  double* target_vec, CSCMatrix<double>*  mat)
{
  // Sanity checks.
  if (target_vec == nullptr) throw Exceptions::NullException(3);
  // Define temporary projection weak form.
  WeakFormSharedPtr<double> proj_wf(new WeakForm<double>(1));
  proj_wf->set_ext(source_meshfn);

  ProjectionLumpedMatrixFormVol* matrix_form = new ProjectionLumpedMatrixFormVol(0, 0);
  ProjectionLumpedVectorFormVol* vector_form = new ProjectionLumpedVectorFormVol(0);
  //ProjectionLumpedVectorFormVol* vector_form = new_ ProjectionLumpedVectorFormVol(0, source_meshfn);
  proj_wf->add_matrix_form(matrix_form);
  proj_wf->add_vector_form(vector_form);
  // Call main function.
  project_internal(space, proj_wf, target_vec, mat);
}




