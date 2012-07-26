
//Mass lumping
EpetraMatrix<double>* massLumping(EpetraMatrix<double>* mass_matrix){
	 EpetraMatrix<double>* lumped_matrix = new EpetraMatrix<double>;   //M_L
	 int size = mass_matrix->get_size();
	lumped_matrix->prealloc(size);
		int i,j;
	 for(i = 0; i<size; i++)   
			lumped_matrix->pre_add_ij(i, i);

	lumped_matrix->alloc();
	
	double mass;
	 int n_entries;
	double* vals = new double[size];
 int* idxs = new int[size];
 
 
 	 for(i = 0; i<size; i++) {
 			mass_matrix->extract_row_copy(i, size, n_entries, vals, idxs);
				for(int indx =0;indx<n_entries;indx++){
						if((mass=vals[indx])!=0) {
								j = idxs[indx];
								if(j<=i){
										lumped_matrix->add(i, i, mass);
										if(j!=i)	lumped_matrix->add(j, j, mass);
								}
						}							
			}
	}
	 
	 lumped_matrix->finish();
	 
	 delete [] vals;
	 delete [] idxs;
	return lumped_matrix;
}



