//c_ij_x
class ConvectionOperator_x : public WeakForm<double>
{
public:
//Constructor
  ConvectionOperator_x(int num_of_equations = 4):  WeakForm<double>(num_of_equations){
    add_matrix_form(new Convection_1(0));
    add_matrix_form(new Convection_1(1));
    add_matrix_form(new Convection_1(2));
    add_matrix_form(new Convection_1(3));


	};
//Destructor
	~ConvectionOperator_x()	{
		for(int i=0; i<this->mfvol.size();i++){
			delete get_mfvol()[i];
		}
		WeakForm<double>::delete_all();
	}
protected:
 class Convection_1 : public MatrixFormVol<double>
  {
  public:
     Convection_1(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
	  Scalar result = Scalar(0);
	  for (int i = 0; i < n; i++)
			result += wt[i] * u->val[i]*v->dx[i];
		
	  return result;

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone() { return new  Convection_1(this->i); }
  };

};
//c_ij_y
class ConvectionOperator_y : public WeakForm<double>
{
public:
  ConvectionOperator_y(int num_of_equations = 4):  WeakForm<double>(num_of_equations){
    add_matrix_form(new Convection_2(0));
    add_matrix_form(new Convection_2(1));
    add_matrix_form(new Convection_2(2));
    add_matrix_form(new Convection_2(3));


	};
	~ConvectionOperator_y()	{
		for(int i=0; i<this->mfvol.size();i++){
			delete get_mfvol()[i];
		}
		WeakForm<double>::delete_all();
	}
protected:
 class Convection_2 : public MatrixFormVol<double>
  {
  public:
     Convection_2(int i) : MatrixFormVol<double>(i, i)  {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
	  Scalar result = Scalar(0);
	  for (int i = 0; i < n; i++)
			result += wt[i] * u->val[i]*v->dy[i];;
		
	  return result;

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
        return Ord(20);
    }

    MatrixFormVol<double>* clone() { return new  Convection_2(this->i); }

  };

};


void calculate_D(double rho_i, double rho_v_x_i, double rho_v_y_i, double rho_energy_i, double e_1, double e_2, 
					double rho_j, double rho_v_x_j, double rho_v_y_j, double rho_energy_j, double kappa, double D[4][4]){

	double v_x_mean = (rho_v_x_i/std::sqrt(rho_i) + rho_v_x_j/std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double v_y_mean = (rho_v_y_i/std::sqrt(rho_i) + rho_v_y_j/std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double H_mean = (QuantityCalculator::enthalpy(rho_i, rho_v_x_i,rho_v_y_i, rho_energy_i, kappa)*std::sqrt(rho_i) + QuantityCalculator::enthalpy(rho_j, rho_v_x_j,rho_v_y_j, rho_energy_j, kappa)*std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
double e_norm = std::sqrt(e_1*e_1+e_2*e_2 );
	double v_n = (v_x_mean*e_1+v_y_mean*e_2)/e_norm;
	double q = 0.5*std::sqrt(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
	double b_2 = (kappa-1)/(c_mean*c_mean);
	double b_1 = b_2*q;

	//double lambda[4] ={fabs(v_n- c_mean),fabs(v_n),fabs(v_n + c_mean),fabs(v_n)};

double lambda_1[4] ={fabs(v_x_mean- c_mean),fabs(v_x_mean),fabs(v_x_mean + c_mean),fabs(v_x_mean)};
double lambda_2[4] ={fabs(v_y_mean- c_mean),fabs(v_y_mean),fabs(v_y_mean + c_mean),fabs(v_y_mean)};

//Eintraege siehe Diss Moeller Appendix C
/*	double R[4][4] = { 1, 1, 1, 0, 
										v_x_mean-c_mean*e_1, v_x_mean, v_x_mean+c_mean*e_1, e_2,
										v_y_mean-c_mean*e_2, v_y_mean, v_y_mean+c_mean*e_2, -e_1,
										H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*e_2-v_y_mean*e_1};

	double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-e_1/c_mean), 0.5*(-b_2*v_y_mean-e_2/c_mean), 0.5*b_2,
										1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
										0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+e_1/c_mean), 0.5*(-b_2*v_y_mean+e_2/c_mean), 0.5*b_2,
										e_1*v_y_mean-e_2*v_x_mean, e_2,  -e_1, 0};
*/


double R_1[4][4]= { 1, 1, 1, 0, 
										v_x_mean-c_mean, v_x_mean, v_x_mean+c_mean, 0,
										v_y_mean, v_y_mean, v_y_mean, -1,
										H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, -v_y_mean};

	double R_2[4][4] = { 1, 1, 1, 0, 
										v_x_mean, v_x_mean, v_x_mean, 1,
										v_y_mean-c_mean, v_y_mean, v_y_mean+c_mean, 0,
										H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean};

	double L_1[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-1/c_mean), 0.5*(-b_2*v_y_mean), 0.5*b_2,
										1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
										0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+1/c_mean), 0.5*(-b_2*v_y_mean), 0.5*b_2,
										v_y_mean, 0,  -1, 0};

	double L_2[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean), 0.5*(-b_2*v_y_mean-1/c_mean), 0.5*b_2,
										1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
										0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean), 0.5*(-b_2*v_y_mean+1/c_mean), 0.5*b_2,
										-1*v_x_mean, 1,  0, 0};

double A_1[4][4];
double A_2[4][4];

	//lambda*L
	for(int i =0;i<4;i++)
			for(int j=0;j<4;j++){
					L_1[i][j] *= lambda_1[i];
					L_2[i][j] *= lambda_2[i];
					A_1[i][j]= 0;
					A_2[i][j]= 0;
			}

//A=RL
for(int k=0;k<4;k++)
	for(int i =0;i<4;i++)
			for(int j=0;j<4;j++){
				A_1[i][k] +=R_1[i][j]*L_1[j][k];
				A_2[i][k] +=R_2[i][j]*L_2[j][k];
	}

	for(int i =0;i<4;i++)
			for(int j=0;j<4;j++)
					D[i][j] = fabs(e_1)*A_1[i][j] + fabs(e_2)*A_2[i][j];
			

}



//artificial Diffusion
UMFPackMatrix<double>* artificialDiffusion(double kappa,double* coeff,Space<double>* space_rho,Space<double>* space_rho_v_x, Space<double>* space_rho_v_y,Space<double>* space_e,int dof_rho, int dof_vel_x, int dof_vel_y, int dof_energy,UMFPackMatrix<double>* matrix_K){

 	ConvectionOperator_x conv_1;
	ConvectionOperator_y conv_2;

	UMFPackMatrix<double>* c_matrix_1 = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double>* c_matrix_2 = new UMFPackMatrix<double> ; 

  DiscreteProblem<double> dp_1(&conv_1, Hermes::vector<const Space<double>*>(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  DiscreteProblem<double> dp_2(&conv_2, Hermes::vector<const Space<double>*>(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  dp_1.assemble(c_matrix_1,NULL,true);
	dp_2.assemble(c_matrix_2);


		int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(space_rho, space_rho_v_x, space_rho_v_y, space_e));

		double* coeff_rho = new double[dof_rho];
		double* coeff_vel_x = new double[dof_rho];
		double* coeff_vel_y = new double[dof_rho];
		double* coeff_energy = new double[dof_rho];


for(int i =0;i<dof_rho; i++){
		coeff_rho[i]=coeff[i];
}
	int k = 0;
for(int i = dof_rho; i<(dof_rho+dof_vel_x);i++){
		coeff_vel_x[k] =coeff[i];
		k++;
}
k = 0;
for(int i = (dof_rho+dof_vel_x); i<(dof_rho+dof_vel_x+dof_vel_y);i++){
		coeff_vel_y[k] =coeff[i];
		k++;
}
k = 0;
for(int i = (dof_rho+dof_vel_x+dof_vel_y); i<ndof;i++){
		coeff_energy[k] =coeff[i];
		k++;
}


	 int size = c_matrix_1->get_size();
	 int nnz = c_matrix_1->get_nnz();
	 UMFPackMatrix<double>* diffusion = new UMFPackMatrix<double>;  
	diffusion->create(size, nnz, c_matrix_1->get_Ap(), c_matrix_1->get_Ai(),c_matrix_1->get_Ax());
//double D[4][4]; 
//diffusion->create(matrix_K->get_size(), matrix_K->get_nnz(), matrix_K->get_Ap(), matrix_K->get_Ai(),matrix_K->get_Ax());
	diffusion->zero();  //matrix = 0

	if(size!=ndof) info("dof != size");
	if(dof_rho!=dof_vel_x)info("rho_dof != dof_vel_x");
	if(dof_rho!=dof_vel_y)info("rho_dof != dof_vel_y");
	if(dof_rho!=dof_energy)info("rho_dof != dof_energy");

	double* Ax_1 = c_matrix_1->get_Ax();
	int* Ai_1 = c_matrix_1->get_Ai();
	int* Ap_1 = c_matrix_1->get_Ap();
	double* Ax_2 = c_matrix_2->get_Ax();
	int* Ai_2 = c_matrix_2->get_Ai();
	int* Ap_2 = c_matrix_2->get_Ap();
	double c_ij_x,c_ij_y, c_ji_x,c_ji_y, d_ij, d_ji, e_ij, e_ji, c_j, c_i, abs_c_ij, abs_c_ji;


	for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
			for(int indx = Ap_1[j]; indx<Ap_1[j+1];indx++){	
				int i = Ai_1[indx];	
				if(j<i){
					if(i>dof_rho) info("i groesser als dof_rho! seltsam!!!!");
					c_ij_x = Ax_1[indx];
					c_ij_y = Ax_2[indx];
					c_ji_x = c_matrix_1->get(j,i);
					c_ji_y = c_matrix_2->get(j,i);  
				c_i = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[i], coeff_vel_x[i], coeff_vel_y[i], coeff_energy[i], kappa)/coeff_rho[i]);
				c_j = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[j], coeff_vel_x[j], coeff_vel_y[j], coeff_energy[j], kappa)/coeff_rho[j]);

//Book
					e_ij = 0.5*Hermes::sqrt((c_ji_x-c_ij_x)*(c_ji_x-c_ij_x)+(c_ji_y-c_ij_y)*(c_ji_y-c_ij_y));
					e_ji = 0.5*Hermes::sqrt((c_ij_x-c_ji_x)*(c_ij_x-c_ji_x)+(c_ij_y-c_ji_y)*(c_ij_y-c_ji_y));
					d_ij = fabs( (c_ji_x-c_ij_x)*coeff_vel_x[j]/(2.*coeff_rho[j])+ (c_ji_y-c_ij_y)*coeff_vel_y[j]/(2.*coeff_rho[j]))+ e_ij*c_j;
					d_ji = fabs( (c_ij_x-c_ji_x)*coeff_vel_x[i]/(2.*coeff_rho[i])+ (c_ij_y-c_ji_y)*coeff_vel_y[i]/(2.*coeff_rho[i]))+ e_ji*c_i;

//paper_failsafe
/*
abs_c_ij = Hermes::sqrt(c_ij_x*c_ij_x + c_ij_y*c_ij_y);
abs_c_ji = Hermes::sqrt(c_ji_x*c_ji_x + c_ji_y*c_ji_y);
d_ij = fabs(c_ij_x*coeff_vel_x[j]/coeff_rho[j] + c_ij_y*coeff_vel_y[j]/coeff_rho[j])+ abs_c_ij*c_j;
d_ji = fabs(c_ji_x*coeff_vel_x[i]/coeff_rho[i] + c_ji_y*coeff_vel_y[i]/coeff_rho[i])+ abs_c_ji*c_i;
*/
					for(int k = 0;k<4;k++){
							int next = k*dof_rho;
						if((i+next>size)||(j+next>size)) info("groesser als size");
						if(d_ij>d_ji){
							diffusion->add(i+next,j+next,d_ij);
							diffusion->add(j+next,i+next,d_ij);	
							diffusion->add(j+next,j+next,-d_ij);
							diffusion->add(i+next,i+next,-d_ij);	
						}else{
							diffusion->add(i+next,j+next,d_ji);
							diffusion->add(j+next,i+next,d_ji);	
							diffusion->add(j+next,j+next,-d_ji);
							diffusion->add(i+next,i+next,-d_ji);
					 }
					}

//Gurris, Moeller, Book (43)
//Diffusion-matrix auf matrix_K setzen!
/*
calculate_D(coeff_rho[i], coeff_vel_x[i] , coeff_vel_y[i],  coeff_energy[i] , 0.5*(c_ij_x-c_ji_x), 0.5*(c_ij_y-c_ji_y), 
					coeff_rho[j], coeff_vel_x[j] , coeff_vel_y[j],  coeff_energy[j], kappa, D);

for(int k = 0;k<4;k++){
		for(int l=k; l<4;l++){
							int next_k = k*dof_rho;
							int next_l = l*dof_rho;
						if((i+next_k>size)||(j+next_l>size)) info("groesser als size");			
						if(D[k][l]!=0.0){
							diffusion->add(i+next_k,j+next_l,D[k][l]);
							diffusion->add(j+next_l,i+next_k,D[k][l]);	
							diffusion->add(j+next_l,j+next_l,-D[k][l]);
							diffusion->add(i+next_k,i+next_k,-D[k][l]);
						}					 
					}
}
*/

				}
	  }
	}

delete [] coeff_rho;
delete [] coeff_vel_x;
delete [] coeff_vel_y;
delete [] coeff_energy;
delete c_matrix_1;
delete c_matrix_2;


	return diffusion;

}

