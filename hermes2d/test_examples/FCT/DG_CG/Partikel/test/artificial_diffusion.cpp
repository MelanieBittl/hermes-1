//c_ij_x
class ConvectionOperator_x : public WeakForm<double>
{
public:
//Constructor
  ConvectionOperator_x(int num_of_equations = 8):  WeakForm<double>(num_of_equations){
for(int i = 0; i<8; i++)
    add_matrix_form(new Convection_1(i));



	};

	
		WeakForm<double>* clone() const
    {
      const_cast<ConvectionOperator_x*>(this)->warned_nonOverride = false;
      return new ConvectionOperator_x(*this);
    }
protected:
 class Convection_1 : public MatrixFormVol<double>
  {
  public:
     Convection_1(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
	  Scalar result = Scalar(0);
	  for (int i = 0; i < n; i++)
			result += wt[i] * u->val[i]*v->dx[i];
		
	  return result;

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone() const { return new  Convection_1(this->i); }
  };

};
//c_ij_y
class ConvectionOperator_y : public WeakForm<double>
{
public:
  ConvectionOperator_y(int num_of_equations = 8):  WeakForm<double>(num_of_equations){
for(int i = 0; i<8; i++)
    add_matrix_form(new Convection_2(i));
	};
	
	WeakForm<double>* clone() const
    {
      const_cast<ConvectionOperator_y*>(this)->warned_nonOverride = false;
      return new ConvectionOperator_y(*this);
    }
	
protected:
 class Convection_2 : public MatrixFormVol<double>
  {
  public:
     Convection_2(int i) : MatrixFormVol<double>(i, i)  {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
	  Scalar result = Scalar(0);
	  for (int i = 0; i < n; i++)
			result += wt[i] * u->val[i]*v->dy[i];;
		
	  return result;

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, Func<double>  **ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>  **ext) const 
    {
        return Ord(20);
    }

    MatrixFormVol<double>* clone() const { return new  Convection_2(this->i); }

  };

};


void calculate_D(double rho_i, double rho_v_x_i, double rho_v_y_i, double rho_energy_i, double e_1, double e_2, 
					double rho_j, double rho_v_x_j, double rho_v_y_j, double rho_energy_j, double gamma, double D[4][4]){

	double v_x_mean = (rho_v_x_i/std::sqrt(rho_i) + rho_v_x_j/std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double v_y_mean = (rho_v_y_i/std::sqrt(rho_i) + rho_v_y_j/std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double H_mean = (QuantityCalculator::enthalpy(rho_i, rho_v_x_i,rho_v_y_i, rho_energy_i, gamma)*std::sqrt(rho_i) + QuantityCalculator::enthalpy(rho_j, rho_v_x_j,rho_v_y_j, rho_energy_j, gamma)*std::sqrt(rho_j))/(std::sqrt(rho_i) +std::sqrt(rho_j));
	double c_mean = std::sqrt((gamma-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
double e_norm = std::sqrt(e_1*e_1+e_2*e_2 );
	double v_n = (v_x_mean*e_1+v_y_mean*e_2)/e_norm;
	double q = 0.5*std::sqrt(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
	double b_2 = (gamma-1)/(c_mean*c_mean);
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
CSCMatrix<double>* artificialDiffusion(double gamma,double* coeff,Hermes::vector<SpaceSharedPtr<double> >  spaces,CSCMatrix<double>* matrix_K)
{

 	ConvectionOperator_x conv_1;
	ConvectionOperator_y conv_2;

	CSCMatrix<double>* c_matrix_1 = new CSCMatrix<double> ; 
	CSCMatrix<double>* c_matrix_2 = new CSCMatrix<double> ; 

  DiscreteProblem<double> dp_1(&conv_1, spaces);
  DiscreteProblem<double> dp_2(&conv_2, spaces);
  dp_1.assemble(c_matrix_1);
	dp_2.assemble(c_matrix_2);


int dof_rho = spaces[0]->get_num_dofs(); Space<double>::assign_dofs(spaces);
int dof_vel_x =dof_rho; int dof_vel_y =dof_rho; int dof_energy =dof_rho;

		int ndof = Space<double>::get_num_dofs(spaces);
		int ndof_g = ndof/2.;

		double* coeff_rho_g = new double[dof_rho];
		double* coeff_vel_x_g = new double[dof_rho];
		double* coeff_vel_y_g = new double[dof_rho];
		double* coeff_e_g = new double[dof_rho];

		double* coeff_rho_p = new double[dof_rho];
		double* coeff_vel_x_p = new double[dof_rho];
		double* coeff_vel_y_p = new double[dof_rho];
		double* coeff_e_p = new double[dof_rho];


for(int i =0;i<dof_rho; i++){
		coeff_rho_g[i]=coeff[i]; coeff_rho_p[i]=coeff[i+ndof_g];
}
	int k = 0;
for(int i = dof_rho; i<(dof_rho+dof_vel_x);i++){
		coeff_vel_x_g[k] =coeff[i]; coeff_vel_x_p[k] =coeff[i+ndof_g];
		k++;
}
k = 0;
for(int i = (dof_rho+dof_vel_x); i<(dof_rho+dof_vel_x+dof_vel_y);i++){
		coeff_vel_y_g[k] =coeff[i]; coeff_vel_y_p[k] =coeff[i+ndof_g];
		k++;
}
k = 0;
for(int i = (dof_rho+dof_vel_x+dof_vel_y); i<ndof_g;i++){
		coeff_e_g[k] =coeff[i]; coeff_e_p[k] =coeff[i+ndof_g];
		k++;
}


	/* int size = c_matrix_1->get_size();
	 int nnz = c_matrix_1->get_nnz();
	 CSCMatrix<double>* diffusion = new CSCMatrix<double>;  
	diffusion->create(size, nnz, c_matrix_1->get_Ap(), c_matrix_1->get_Ai(),c_matrix_1->get_Ax());
	diffusion->zero();  //matrix = 0*/

 int size = matrix_K->get_size();
	 int nnz = matrix_K->get_nnz();
	 CSCMatrix<double>* diffusion = new CSCMatrix<double>;  
	diffusion->create(size, nnz, matrix_K->get_Ap(), matrix_K->get_Ai(), matrix_K->get_Ax());
	diffusion->zero(); 

	if(size!=ndof) printf("dof != size");
	if(dof_rho!=dof_vel_x)printf("rho_dof != dof_vel_x");
	if(dof_rho!=dof_vel_y)printf("rho_dof != dof_vel_y");
	if(dof_rho!=dof_energy)printf("rho_dof != dof_energy");

	double* Ax_1 = c_matrix_1->get_Ax();
	int* Ai_1 = c_matrix_1->get_Ai();
	int* Ap_1 = c_matrix_1->get_Ap();
	double* Ax_2 = c_matrix_2->get_Ax();
	int* Ai_2 = c_matrix_2->get_Ai();
	int* Ap_2 = c_matrix_2->get_Ap();
	double c_ij_x,c_ij_y, c_ji_x,c_ji_y, d_ij, d_ji, e_ij, e_ji, c_j, c_i, abs_c_ij, abs_c_ji;
	double k_ij,k_ji, d_ij_p, d_ij_g;


	for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
			for(int indx = Ap_1[j]; indx<Ap_1[j+1];indx++){	
				int i = Ai_1[indx];	
				if(j<i){
					if(i>dof_rho) printf("i groesser als dof_rho! seltsam!!!!");

//Gasphase
					c_ij_x = Ax_1[indx];
					c_ij_y = Ax_2[indx];
					c_ji_x = c_matrix_1->get(j,i);
					c_ji_y = c_matrix_2->get(j,i); 
				c_i = Hermes::sqrt(gamma*QuantityCalculator::calc_pressure(coeff_rho_g[i], coeff_vel_x_g[i], coeff_vel_y_g[i], coeff_e_g[i], gamma)/coeff_rho_g[i]);
				c_j = Hermes::sqrt(gamma*QuantityCalculator::calc_pressure(coeff_rho_g[j], coeff_vel_x_g[j], coeff_vel_y_g[j], coeff_e_g[j], gamma)/coeff_rho_g[j]);

//Book
					e_ij = 0.5*Hermes::sqrt((c_ji_x-c_ij_x)*(c_ji_x-c_ij_x)+(c_ji_y-c_ij_y)*(c_ji_y-c_ij_y));
					e_ji = 0.5*Hermes::sqrt((c_ij_x-c_ji_x)*(c_ij_x-c_ji_x)+(c_ij_y-c_ji_y)*(c_ij_y-c_ji_y));
					d_ij = fabs( (c_ji_x-c_ij_x)*coeff_vel_x_g[j]/(2.*coeff_rho_g[j])+ (c_ji_y-c_ij_y)*coeff_vel_y_g[j]/(2.*coeff_rho_g[j]))+ e_ij*c_j;
					d_ji = fabs( (c_ij_x-c_ji_x)*coeff_vel_x_g[i]/(2.*coeff_rho_g[i])+ (c_ij_y-c_ji_y)*coeff_vel_y_g[i]/(2.*coeff_rho_g[i]))+ e_ji*c_i;

					d_ij_g = std::max(d_ij,d_ji);



					for(int k = 0;k<4;k++){
							int next = k*dof_rho;
						if((i+next>ndof_g)||(j+next>ndof_g)) printf("gas phase diffusion exceeds range");						
							diffusion->add(i+next,j+next,d_ij_g);
							diffusion->add(j+next,i+next,d_ij_g);	
							diffusion->add(j+next,j+next,-d_ij_g);
							diffusion->add(i+next,i+next,-d_ij_g);			
					}

//Particulate-Phase

				k_ij = std::fabs(c_ij_x *coeff_vel_x_p[j]/coeff_rho_p[j] + c_ij_y *coeff_vel_y_p[j]/coeff_rho_p[j]);
				k_ji = std::fabs(c_ji_x *coeff_vel_x_p[i]/coeff_rho_p[i] + c_ji_y *coeff_vel_y_p[i]/coeff_rho_p[i]);
				d_ij_p = std::max(d_ij,d_ji);
					for(int k = 0;k<4;k++){
							int next = k*dof_rho+ndof_g;
						if((i+next>size)||(j+next>size)) printf("particle phase diffusion exceeds range");						
							diffusion->add(i+next,j+next,d_ij_p);
							diffusion->add(j+next,i+next,d_ij_p);	
							diffusion->add(j+next,j+next,-d_ij_p);
							diffusion->add(i+next,i+next,-d_ij_p);			
					}

				}
	  }
	}

delete c_matrix_1; 
delete c_matrix_2;


	delete [] coeff_rho_g;
	delete [] coeff_vel_x_g;
	delete [] coeff_vel_y_g;
	delete [] coeff_e_g; 

	return diffusion;

}

