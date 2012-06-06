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



//artificial Diffusion
//UMFPackMatrix<double>* artificialDiffusion(double kappa,Solution<double>* prev_rho, Solution<double>* prev_rho_vel_x, Solution<double>* prev_rho_vel_y, Solution<double>* prev_energy,H1Space<double>* space_rho,H1Space<double>* space_rho_v_x, H1Space<double>* space_rho_v_y,H1Space<double>* space_e){

UMFPackMatrix<double>* artificialDiffusion(double kappa,double* coeff,H1Space<double>* space_rho,H1Space<double>* space_rho_v_x, H1Space<double>* space_rho_v_y,H1Space<double>* space_e){

 	ConvectionOperator_x conv_1;
	ConvectionOperator_y conv_2;

	UMFPackMatrix<double>* c_matrix_1 = new UMFPackMatrix<double> ; 
	UMFPackMatrix<double>* c_matrix_2 = new UMFPackMatrix<double> ; 

  DiscreteProblem<double> dp_1(&conv_1, Hermes::vector<const Space<double>*>(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  DiscreteProblem<double> dp_2(&conv_2, Hermes::vector<const Space<double>*>(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  dp_1.assemble(c_matrix_1,NULL,true);
	dp_2.assemble(c_matrix_2);

		int dof_rho = space_rho->get_num_dofs();
		int dof_vel_x = space_rho_v_x->get_num_dofs();
		int dof_vel_y = space_rho_v_y->get_num_dofs();
		int dof_energy = space_e->get_num_dofs();
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

/*Lumped_Projection::project_lumped(space_rho, prev_rho, coeff_rho, matrix_solver);
		Lumped_Projection::project_lumped(space_rho_v_x, prev_rho_vel_x, coeff_vel_x, matrix_solver);
		Lumped_Projection::project_lumped(space_rho_v_y, prev_rho_vel_y, coeff_vel_y, matrix_solver);
		Lumped_Projection::project_lumped(space_e, prev_energy, coeff_energy, matrix_solver);

		/*		OGProjection<double>::project_global(space_rho, prev_rho, coeff_rho, matrix_solver, HERMES_L2_NORM);
		OGProjection<double>::project_global(space_rho_v_x, prev_rho_vel_x, coeff_vel_x, matrix_solver, HERMES_L2_NORM);
		OGProjection<double>::project_global(space_rho_v_y, prev_rho_vel_y, coeff_vel_y, matrix_solver, HERMES_L2_NORM);
		OGProjection<double>::project_global(space_e, prev_energy, coeff_energy, matrix_solver, HERMES_L2_NORM);*/


	 int size = c_matrix_1->get_size();
	 int nnz = c_matrix_1->get_nnz();
	 UMFPackMatrix<double>* diffusion = new UMFPackMatrix<double>;  
	diffusion->create(size, nnz, c_matrix_1->get_Ap(), c_matrix_1->get_Ai(),c_matrix_1->get_Ax());
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
	double c_ij_x,c_ij_y, c_ji_x,c_ji_y, d_ij, d_ji, e_ij, e_ji, c_j, c_i;

//fabs(x)

		for(int j = 0; j<dof_rho; j++){ //Spalten durchlaufen
				for(int indx = Ap_1[j]; indx<Ap_1[j+1];indx++){	
					int i = Ai_1[indx];	
					if(j<i){
						if(i>dof_rho) info("i groesser als dof_rho! seltsam!!!!");
						c_ij_x = Ax_1[indx];
						c_ij_y = Ax_2[indx];
						c_ji_x = c_matrix_1->get(j,i);
						c_ji_y = c_matrix_2->get(j,i);  
						e_ij = 0.5*Hermes::sqrt((c_ji_x-c_ij_x)*(c_ji_x-c_ij_x)+(c_ji_y-c_ij_y)*(c_ji_y-c_ij_y));
						e_ji = 0.5*Hermes::sqrt((c_ij_x-c_ji_x)*(c_ij_x-c_ji_x)+(c_ij_y-c_ji_y)*(c_ij_y-c_ji_y));
					c_i = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[i], coeff_vel_x[i], coeff_vel_y[i], coeff_energy[i], kappa)/coeff_rho[i]);
				c_j = Hermes::sqrt(kappa*QuantityCalculator::calc_pressure(coeff_rho[j], coeff_vel_x[j], coeff_vel_y[j], coeff_energy[j], kappa)/coeff_rho[j]);
					d_ij= fabs( (c_ji_x-c_ij_x)*coeff_vel_x[j]/(2.*coeff_rho[j])+ (c_ji_y-c_ij_y)*coeff_vel_y[j]/(2.*coeff_rho[j]))+ e_ij*c_j;
		d_ji= fabs( (c_ij_x-c_ji_x)*coeff_vel_x[i]/(2.*coeff_rho[i])+ (c_ij_y-c_ji_y)*coeff_vel_y[i]/(2.*coeff_rho[i]))+ e_ji*c_i;

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
				}
	  }
	}


	return diffusion;

}

