#include "lumped_projection.h"
#include "discrete_problem.h"




////template<typename double>
void Lumped_Projection::project_internal( SpaceSharedPtr<double> space, WeakForm<double>* wf, double* target_vec,
                               MatrixSolverType matrix_solver, CSCMatrix<double>*  mat)
{
            // Sanity check.
      if(wf == NULL)
        throw Hermes::Exceptions::NullException(1);
      if(target_vec == NULL)
        throw Exceptions::NullException(2);
    if(space == NULL)  throw Hermes::Exceptions::Exception("this->space == NULL in project_internal().");


      // Get dimension of the space.
      int ndof = space->get_num_dofs();

if(mat!=NULL) if(mat->get_size()!=ndof) printf("lumped_projection: matrixsize =%i !=  ndof=%i", mat->get_size(),ndof);

  // Initialize DiscreteProblem.
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(wf, space);
      dp->set_linear();
      dp->set_do_not_use_cache();

    	CSCMatrix<double>* matrix = new CSCMatrix<double>;	
  	SimpleVector<double>* rhs = new SimpleVector<double>(ndof);
	double* coeff_vec =NULL; 

	if(mat==NULL) 
	{ 
		CSCMatrix<double>* lumped_matrix = new CSCMatrix<double>;   //M_L 
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
  try
  {
   solver->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
			coeff_vec = solver->get_sln_vector();		

		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
		delete lumped_matrix;

	}else{ 
		dp->assemble(rhs);
		UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(mat,rhs);		
  try
  {
   solver->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	
			coeff_vec = solver->get_sln_vector();			
		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
	}
  
  
  delete matrix;
  delete rhs;
  delete dp;

  
}


  void Lumped_Projection::project_lumped(Hermes::vector<SpaceSharedPtr<double> >  spaces, 
        Hermes::vector<MeshFunctionSharedPtr<double> > source_meshfns,
        double* target_vec, Hermes::MatrixSolverType matrix_solver, CSCMatrix<double>*  mat)
    {
     
      int n = spaces.size();

      // Sanity checks.
      if (n != source_meshfns.size()) throw Exceptions::LengthException(1, 2, n, source_meshfns.size());
      if (target_vec == NULL) throw Exceptions::NullException(3);
  
      int start_index = 0;
      for (int i = 0; i < n; i++) 
			{     
				spaces[i]->assign_dofs();
          project_lumped(spaces[i], source_meshfns[i], target_vec + start_index, matrix_solver, mat);
        start_index += spaces[i]->get_num_dofs();
			}
		Space<double>::assign_dofs(spaces);
    }


void Lumped_Projection::project_lumped( SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double>  source_meshfn,
                             double* target_vec, MatrixSolverType matrix_solver ,CSCMatrix<double>*  mat )
{
			

      // Sanity checks.
      if (target_vec == NULL) throw Exceptions::NullException(3);

      // Define temporary projection weak form.
      WeakForm<double>* proj_wf = new WeakForm<double>(1);
      proj_wf->set_ext(source_meshfn);
      // Add Jacobian.
      proj_wf->add_matrix_form(new MatrixDefaultNormFormVol<double>(0, 0, HERMES_L2_NORM));
      // Add Residual.
      proj_wf->add_vector_form(new VectorDefaultNormFormVol<double>(0, HERMES_L2_NORM));

      // Call main function.
      project_internal(space, proj_wf, target_vec, matrix_solver, mat);

      // Clean up.
      delete proj_wf;


}


