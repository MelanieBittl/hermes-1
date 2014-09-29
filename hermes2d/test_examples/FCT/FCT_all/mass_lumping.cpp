//Mass lumping
CSCMatrix<double>* massLumping(CSCMatrix<double>* mass_matrix){
	 CSCMatrix<double>* lumped_matrix = new CSCMatrix<double>;   //M_L
	 int size = mass_matrix->get_size();
	 double diag[size];
	 int nnz = mass_matrix->get_nnz();
	 int row[size]; 
	int col[size+1];
int i;
	 for(i = 0; i<size; i++){    
	    diag[i] = 0;
	    row[i]= i;
	    col[i]=i;
	 }
	col[size]=size;// = Anzahl der Nichtnulleintraege (bezieht sich auf theoretisch als naechstes kommenden Eintrag)
	 for(i = 0; i<nnz; i++){    
	    diag[mass_matrix->get_Ai()[i]] += mass_matrix->get_Ax()[i]; 
	 }
	 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
	return lumped_matrix;
}


