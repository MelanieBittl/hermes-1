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
  ConvectionOperator_y(int num_of_equations = 4):  WeakForm<double>(num_of_equations){
    add_matrix_form(new Convection_2(0));
    add_matrix_form(new Convection_2(1));
    add_matrix_form(new Convection_2(2));
    add_matrix_form(new Convection_2(3));


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
			result += wt[i] * u->val[i]*v->dy[i];
		
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



//test
class ConvectionOperator_test : public WeakForm<double>
{
public:
  ConvectionOperator_test(int num_of_equations = 4):  WeakForm<double>(num_of_equations){
    add_matrix_form(new Convection_0(0));
    add_matrix_form(new Convection_0(1));
    add_matrix_form(new Convection_0(2));
    add_matrix_form(new Convection_0(3));


	};
	
	WeakForm<double>* clone() const
    {
      const_cast<ConvectionOperator_test*>(this)->warned_nonOverride = false;
      return new ConvectionOperator_test(*this);
    }
	
protected:
 class Convection_0 : public MatrixFormVol<double>
  {
  public:
     Convection_0(int i) : MatrixFormVol<double>(i, i)  {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, Func<Scalar>  **ext) const 
    {
	  Scalar result_x = Scalar(0);
	  Scalar result_y = Scalar(0);
	  for (int i = 0; i < n; i++)
			result_x += wt[i] * u->val[i]*v->dx[i];
	  for (int i = 0; i < n; i++)
			result_y += wt[i] * u->val[i]*v->dy[i];
		
if(result_x!=0.) return result_x;
else if(result_y!=0.) return result_y;
else	  return 0.;

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

    MatrixFormVol<double>* clone() const { return new  Convection_0(this->i); }

  };

};


/////////////////////////////////////////////////////////////////////////////////
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


//-----------------------------------------------
//--------------artificial Diffusion--------------
//-----------------------------------------------
CSCMatrix<double>* artificialDiffusion(double kappa,double* coeff,Hermes::vector<SpaceSharedPtr<double> >  spaces,CSCMatrix<double>* mat_KS, bool* fct = NULL){

 	ConvectionOperator_x conv_1;
	ConvectionOperator_y conv_2;


	CSCMatrix<double>* c_matrix_1 = new CSCMatrix<double> ; 
	CSCMatrix<double>* c_matrix_2 = new CSCMatrix<double> ; 

  DiscreteProblem<double> dp_1(&conv_1, spaces);
  DiscreteProblem<double> dp_2(&conv_2, spaces);
  dp_1.assemble(c_matrix_1);
	dp_2.assemble(c_matrix_2);


		int ndof = Space<double>::get_num_dofs(spaces);
int dof_rho = spaces[0]->get_num_dofs();
 int dof_vel_x= spaces[1]->get_num_dofs();
int dof_vel_y= spaces[2]->get_num_dofs(); 
int dof_energy= spaces[3]->get_num_dofs();


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


	/* int size = c_matrix_1->get_size();
	 int nnz = c_matrix_1->get_nnz();*/
	 CSCMatrix<double>* diffusion = new CSCMatrix<double>; 
	//ConvectionOperator_test conv_0;
 // DiscreteProblem<double> dp_3(&conv_0, spaces);
	//dp_3.assemble(diffusion);
	//diffusion->create(size, nnz, c_matrix_1->get_Ap(), c_matrix_1->get_Ai(),c_matrix_1->get_Ax());
	 int size = mat_KS->get_size();
	 int nnz = mat_KS->get_nnz();
diffusion->create_pattern(size, nnz, mat_KS->get_Ap(), mat_KS->get_Ai()); //matrix = 0


	if(size!=ndof)   	Hermes::Mixins::Loggable::Static::info("dof != size");
	if(dof_rho!=dof_vel_x)  	Hermes::Mixins::Loggable::Static::info("rho_dof != dof_vel_x");
	if(dof_rho!=dof_vel_y)  	Hermes::Mixins::Loggable::Static::info("rho_dof != dof_vel_y");
	if(dof_rho!=dof_energy)  	Hermes::Mixins::Loggable::Static::info("rho_dof != dof_energy");

	double* Ax_1 = c_matrix_1->get_Ax();
	int* Ai_1 = c_matrix_1->get_Ai();
	int* Ap_1 = c_matrix_1->get_Ap();
	double* Ax_2 = c_matrix_2->get_Ax();
	int* Ai_2 = c_matrix_2->get_Ai();
	int* Ap_2 = c_matrix_2->get_Ap();
	double c_ij_x,c_ij_y, c_ji_x,c_ji_y, d_ij, d_ji, e_ij, e_ji, c_j, c_i, abs_c_ij, abs_c_ji;

	double a = 0.;
	for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
		if((fct!=NULL)&&(fct[j]==false)) continue;
			for(int indx = Ap_1[j]; indx<Ap_1[j+1];indx++){	
				int i = Ai_1[indx];
				if((fct!=NULL)&&(fct[i]==false)) continue;	
				if(j<i){
					if(i>dof_rho)   	Hermes::Mixins::Loggable::Static::info("i groesser als dof_rho! seltsam!!!!");
					c_ij_x = Ax_1[indx];
					c_ij_y =	c_matrix_2->get(i,j);  
					c_ji_x = c_matrix_1->get(j,i);
					c_ji_y = c_matrix_2->get(j,i);  
				c_i = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[i], coeff_vel_x[i], coeff_vel_y[i], coeff_energy[i], kappa)/coeff_rho[i]);
				c_j = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[j], coeff_vel_x[j], coeff_vel_y[j], coeff_energy[j], kappa)/coeff_rho[j]);

//Book
					e_ij = 0.5*Hermes::sqrt((c_ji_x-c_ij_x)*(c_ji_x-c_ij_x)+(c_ji_y-c_ij_y)*(c_ji_y-c_ij_y));
					e_ji = 0.5*Hermes::sqrt((c_ij_x-c_ji_x)*(c_ij_x-c_ji_x)+(c_ij_y-c_ji_y)*(c_ij_y-c_ji_y));
					d_ij = fabs( (c_ji_x-c_ij_x)*coeff_vel_x[j]/(2.*coeff_rho[j])+ (c_ji_y-c_ij_y)*coeff_vel_y[j]/(2.*coeff_rho[j]))+ e_ij*c_j;
					d_ji = fabs( (c_ij_x-c_ji_x)*coeff_vel_x[i]/(2.*coeff_rho[i])+ (c_ij_y-c_ji_y)*coeff_vel_y[i]/(2.*coeff_rho[i]))+ e_ji*c_i;

					if(d_ij>d_ji)
						 a = d_ij; 
					else a = d_ji;
					if(a!=0.)
					{
							for(int k = 0;k<4;k++){
									int next = k*dof_rho;
								if((i+next>size)||(j+next>size))   	Hermes::Mixins::Loggable::Static::info("groesser als size");
									diffusion->add(i+next,j+next,a);
									diffusion->add(j+next,i+next,a);
									diffusion->add(j+next,j+next,-a);
									diffusion->add(i+next,i+next,-a);
								}	
					 }				
				}
	  }
	}

delete c_matrix_1; 
delete c_matrix_2;


	delete [] coeff_rho;
	delete [] coeff_vel_x;
	delete [] coeff_vel_y;
	delete [] coeff_energy; 

	return diffusion;

}

