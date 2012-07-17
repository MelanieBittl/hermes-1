#include "lumped_projection.h"



#include "discrete_problem.h"



////template<typename double>
void Lumped_Projection::project_internal( const Space<double>* space, WeakForm<double>* wf, double* target_vec,
                               MatrixSolverType matrix_solver, UMFPackMatrix<double>*  mat)
{

      // Sanity check.
    if(space == NULL) throw Hermes::Exceptions::Exception("this->space == NULL in project_internal().");

      // Get dimension of the space.
      int ndof = space->get_num_dofs();

	if(mat!=NULL) if(mat->get_size()!=ndof) printf("matrix=%i, ndof=%i", mat->get_size(),ndof);

  // Initialize DiscreteProblem.
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(wf, space);
    	UMFPackMatrix<double>* matrix = new UMFPackMatrix<double>;	
  	UMFPackVector<double>* rhs = new UMFPackVector<double>(ndof);
	double* coeff_vec =NULL; 
	if(mat==NULL) { 
		UMFPackMatrix<double>* lumped_matrix = new UMFPackMatrix<double>;   //M_L 
		dp->assemble(matrix, rhs);  
			//Masslumping		 
		 int size = matrix->get_size();
		 double diag[size];
		 int nnz = matrix->get_nnz();
		 int row[size]; 
		int col[size+1];
		 for(int i = 0; i<size; i++){    
		    diag[i] = 0;
		    row[i]= i;
		    col[i]=i;
		 }
		col[size]=size;// letzter Eintrag bezieht sich auf nichts(Ende der Matrix)= Anzahl der Eintraege!
	
		 for(int i = 0; i<nnz; i++){    
		    diag[matrix->get_Ai()[i]] += matrix->get_Ax()[i]; 
		 }

		 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
		UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(lumped_matrix,rhs);		
		if(solver->solve()){ 
			coeff_vec = solver->get_sln_vector();		
		}
	  	else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
		delete lumped_matrix;
	}else{ 
		dp->assemble(rhs);
		UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(mat,rhs);		
		if(solver->solve()) 
			coeff_vec = solver->get_sln_vector();			
	 	 else throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
	}
  
  
  delete matrix;
  delete rhs;
  delete dp;
  
  
}

void Lumped_Projection::project_lumped( const  Space<double>* space, MeshFunction<double>* source_meshfn,
                             double* target_vec, MatrixSolverType matrix_solver ,UMFPackMatrix<double>*  mat )
{
		
      // Sanity checks.
      if (target_vec == NULL) throw Exceptions::NullException(3);

      // Define temporary projection weak form.
      WeakForm<double>* proj_wf = new WeakForm<double>(1);

			ProjectionLumpedMatrixFormVol* matrix_form =	new ProjectionLumpedMatrixFormVol(0, 0);
			ProjectionLumpedVectorFormVol* vector_form = new ProjectionLumpedVectorFormVol(0, source_meshfn);
      // Add Jacobian.
      proj_wf->add_matrix_form(matrix_form);
      // Add Residual.
      proj_wf->add_vector_form(vector_form);

      // Call main function.
      project_internal(space, proj_wf, target_vec, matrix_solver, mat);

      // Clean up.
      delete proj_wf;
			delete vector_form;
			delete matrix_form;

}


