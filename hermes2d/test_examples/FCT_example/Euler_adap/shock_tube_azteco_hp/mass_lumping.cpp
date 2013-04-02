
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

//Mass lumping an den Stellen von fct==true, sonst standard Massmatrix
EpetraMatrix<double>* massLumping(EpetraMatrix<double>* mass_matrix,	DiscreteProblem<double>* dp_mass,bool* fct)
{  	
	int size = mass_matrix->get_size();
	EpetraMatrix<double>* lumped_matrix = new EpetraMatrix<double>;   //M_L
 	dp_mass->assemble(lumped_matrix);
	double a =0.0; int j,n_entries;
	double* vals = new double[size];
	int* idxs = new int[size]; 
	/*	for(int j = 0; j<size; j++){ //Spalten durchlaufen
				if(fct[j]== false) continue;
				for(int i=0; i<size;i++){	
							if(fct[i]== false) continue;
							if(j<i){
								a = mass_matrix->get(i,j);
								if(a!=0.0){
									lumped_matrix->add(i,i,a);    //zur Diagonale hinzufuegen
									lumped_matrix->add(j,j,a);    //zur Diagonale hinzufuegen
									lumped_matrix->add(i,j,-a);	  //i,j Eintrag auf 0 setzen
									lumped_matrix->add(j,i,-a);	 //j,i Eintrag auf 0 setzen
								}	
							}
				}
			}*/
			
	for(int i = 0; i<size; i++) {
		if(fct[i]== false) continue;
		mass_matrix->extract_row_copy(i, size, n_entries, vals, idxs);
		for(int indx =0;indx<n_entries;indx++){
			if((a=vals[indx])!=0) {
				j = idxs[indx];
				if(fct[j]== false) continue;
				if(j<i){
					if(a!=0.0){
						lumped_matrix->add(i,i,a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(j,j,a);    //zur Diagonale hinzufuegen
						lumped_matrix->add(i,j,-a);	  //i,j Eintrag auf 0 setzen
						lumped_matrix->add(j,i,-a);	 //j,i Eintrag auf 0 setzen
					}
				}
			}							
		}
	}
			
	 delete [] vals;
	 delete [] idxs;		
	return lumped_matrix;
}


