 
void calc_z_z_error( Space<double>* ref_space, Solution<double>* u, Solution<double>* R_h_1, Solution<double>* R_h_2, double* elem_error)
{
	int ndof = ref_space->get_num_dofs();
	GradientReconstruction_1* grad_1 = new GradientReconstruction_1(u);
	GradientReconstruction_2* grad_2 = new GradientReconstruction_2(u);

	DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, ref_space);
	DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, ref_space);
	UMFPackMatrix<double> * mass_matrix = new UMFPackMatrix<double> ; 
	UMFPackVector<double> * rhs_1 = new UMFPackVector<double>(ndof);
	UMFPackVector<double> * rhs_2 = new UMFPackVector<double>(ndof);
	dp_1->assemble(mass_matrix,rhs_1); 
	dp_2->assemble(rhs_2);
	UMFPackLinearSolver<double> * solver_1 = new UMFPackLinearSolver<double> (mass_matrix,rhs_1);
	if(solver_1->solve()){ 				
		Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , ref_space, R_h_1);	
	}else error ("Matrix solver failed.\n");
	UMFPackLinearSolver<double> * solver_2 = new UMFPackLinearSolver<double> (mass_matrix,rhs_2);	
	if(solver_2->solve()){ 				
		Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , ref_space, R_h_2);	
	}else error ("Matrix solver failed.\n");



      // order of integral
      const int order = 20;

      // initialize total_error
      double total_error = 0.0;

      // element for the element loop
      Element *e =NULL;

      // refmap for computing Jacobian
			RefMap* rm = new RefMap;
      rm->set_quad_2d(&g_quad_2d_std);

      // loop over elements
			int id;

      for_all_active_elements( e, ref_space->get_mesh() ) {

					elem_error[e->id] =0.;

        // set up the solution quadrature
        u->set_quad_2d(&g_quad_2d_std);
        u->set_active_element(e);
        u->set_quad_order(order);

				
				R_h_1->set_quad_2d(&g_quad_2d_std);
        R_h_1->set_active_element(e);
        R_h_1->set_quad_order(order);
			  R_h_2->set_quad_2d(&g_quad_2d_std);
        R_h_2->set_active_element(e);
        R_h_2->set_quad_order(order);

        // get the constant Jacobian
        rm->set_active_element(e);
        double jac = rm->get_const_jacobian();

        // get the quadrature points
        int np = u->get_quad_2d()->get_num_points(order,HERMES_MODE_QUAD);
        double3 *pt = g_quad_2d_std.get_points(order,HERMES_MODE_QUAD);

        // get the function derivative values
        Func<double>* v = init_fn( u, order );
				Func<double>* R_1 = init_fn( R_h_1, order );
				Func<double>* R_2 = init_fn( R_h_2, order );

        // loop over points and integrate the energy

        for( int j = 0; j < np; ++j ) {
					elem_error[e->id] +=pt[j][2]*jac*( (R_1->val[j]-v->dx[j])*(R_1->val[j]-v->dx[j]) + (R_2->val[j]-v->dy[j])*(R_2->val[j]-v->dy[j]) );          
        }

				total_error += elem_error[e->id];

			v->free_fn();
			R_1->free_fn();
			R_2->free_fn();
				delete v;
				delete R_1;
				delete R_2;
			
      }


		delete rm;

	delete grad_1;
	delete grad_2;
	delete dp_1;
	delete dp_2;
	delete mass_matrix;
	delete rhs_1;
	delete rhs_2;
	delete solver_1;
	delete solver_2;


 }

