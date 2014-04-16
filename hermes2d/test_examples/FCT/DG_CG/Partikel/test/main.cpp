#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "interface.h"
#include "euler_flux.h"
#include "lumped_projection.h"


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 0.000001;
const double T_FINAL = 50000;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// GAMMA.
const double GAMMA = 1.4; 

// Penalty Parameter.
double SIGMA =std::pow(10,8);

//Particle density_particle
const double density_particle = 4000.; 

const double diameter = 2e-5;		
const double c_vg = 743.;
const double c_vp = 1380.;
const double c_pg = 1040;
const double Pr = 0.75;		
const double mu = 2.76*1e-5;


    

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;


void assemble_vector_s(SpaceSharedPtr<double> space,Hermes::vector<MeshFunctionSharedPtr<double> > slns, int dof_total, int dof_gas, int dof_rho, double* s, CSCMatrix<double>* matrix, CSCMatrix<double>* lumped_matrix)
{
	Element* e =NULL;
	AsmList<double> al;
	double x, y; int dof; double F_D_1, F_D_2, Q_T;
matrix->zero();
double kap = c_pg*mu/Pr;
bool visited[dof_gas]; for(int i= 0; i<dof_gas; i++) visited[i] = false; 

for(int i= 0; i<dof_total; i++) s[i] = 0.; 
	
	for_all_active_elements(e, space->get_mesh()){	
		space->get_element_assembly_list(e, &al);
	  	for (unsigned int iv = 0; iv < e->get_nvert(); iv++){   		
		  int index =  space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
			Node* vn = e->vn[iv];
			if (space->get_element_order(e->id) == 0) break;
			if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
				for(unsigned int j = 0; j < al.get_cnt(); j ++){			 
					if((al.get_idx()[j]==index)&&(( dof=al.get_dof()[j])!=-1.0)){ 					
						x = vn->x; 
						y = vn->y;
					if(visited[dof]==false){
						visited[dof]=true;
						double alpha_p  =slns[4]->get_pt_value(x,y,false,e)->val[0]/density_particle ;
						double alpha_g = 1.-alpha_p;

						double rho_g = slns[0]->get_pt_value(x,y,false,e)->val[0]/alpha_g;  
						double rho_v_x_g = slns[1]->get_pt_value(x,y,false,e)->val[0]/alpha_g; 
						double rho_v_y_g = slns[2]->get_pt_value(x,y,false,e)->val[0]/alpha_g; 
						double rho_e_g = slns[3]->get_pt_value(x,y,false,e)->val[0]/alpha_g;
						double v_x_g = rho_v_x_g/rho_g;
						double v_y_g = rho_v_y_g/rho_g;

						double rho_p = density_particle;  
						double rho_v_x_p = slns[5]->get_pt_value(x,y,false,e)->val[0]/alpha_p; 
						double rho_v_y_p = slns[6]->get_pt_value(x,y,false,e)->val[0]/alpha_p; 
						double rho_e_p = slns[7]->get_pt_value(x,y,false,e)->val[0]/alpha_p;
						double v_x_p = rho_v_x_p/rho_p;
						double v_y_p = rho_v_y_p/rho_p;

						double v1_diff = (v_x_g - v_x_p);
						double v2_diff = (v_y_g - v_y_p);
						double v_diff_abs = std::sqrt(v1_diff*v1_diff+ v2_diff*v2_diff);

						double Re = rho_g*diameter*v_diff_abs/mu;
						double C_D = 0.44;
						if(Re<1000) C_D=24./Re*(1.+0.15*std::pow(Re,0.687));
						double Nu = 2.+0.65*std::sqrt(Re)*std::pow(Pr,1./3.);

						double T_g = 1./c_vg*(rho_e_g/rho_g-0.5*(v_x_g*v_x_g+v_y_g*v_y_g));
						double T_p = 1./c_vp*(rho_e_p/rho_p-0.5*(v_x_p*v_x_p+v_y_p*v_y_p));

						F_D_1 = rho_g*alpha_p*C_D/diameter*v1_diff*v_diff_abs*0.75;
						F_D_2 = rho_g*alpha_p*C_D/diameter*v2_diff*v_diff_abs*0.75;
						Q_T = Nu*6.*kap/sqr(diameter)*alpha_p*(T_g-T_p);

						int ind = dof+dof_rho;
						s[ind] = -F_D_1;s[ind+dof_gas] = -s[ind]; ind += dof_rho;
						s[ind] = -F_D_2;s[ind+dof_gas] = -s[ind];ind += dof_rho;
						s[ind] = -(F_D_1*v_x_p+F_D_2*v_y_p)-Q_T; s[ind+dof_gas] = -s[ind];	

						double Q_drag = 0.75*v_diff_abs*C_D/(diameter*alpha_g*density_particle);
						double Q_tem = Nu*6.*kap/(diameter*diameter*density_particle);


						double alpha_rho_g= rho_g*alpha_g; 
						double alpha_rho_v_x_g = rho_v_x_g*alpha_g; 
						double alpha_rho_v_y_g= rho_v_y_g*alpha_g; 
						double alpha_rho_e_g= rho_e_g*alpha_g; 
						double alpha_rho_p = rho_p*alpha_p; 
						double alpha_rho_v_x_p= rho_v_x_p*alpha_p; 
						double alpha_rho_v_y_p = rho_v_y_p*alpha_p; 
						double alpha_rho_e_p = rho_e_p*alpha_p; 

//1. Zeile 0
//2.Zeile
						int ind_i = dof+dof_rho; int ind_j = dof; int ind_i_p = ind_i+dof_gas;
double lumped = lumped_matrix->get(ind_i, ind_i);
						matrix->add(ind_i,ind_j, Q_drag*alpha_rho_v_x_p*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*alpha_rho_v_x_p*lumped); ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_p*lumped); matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_p*lumped);  ind_j = dof+dof_gas;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_v_x_g*lumped);matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_v_x_g*lumped); ind_j += dof_rho;
						matrix->add(ind_i,ind_j, Q_drag*alpha_rho_g*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*alpha_rho_g*lumped); 
//3. Zeile
						ind_i += dof_rho; ind_j = dof; ind_i_p = ind_i+dof_gas;
						matrix->add(ind_i,ind_j, Q_drag*alpha_rho_v_y_p*lumped);matrix->add(ind_i_p,ind_j, -Q_drag*alpha_rho_v_y_p*lumped); ind_j += 2.*dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_p*lumped);matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_p*lumped); ind_j = dof+dof_gas;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_v_y_g*lumped);matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_v_y_g*lumped); ind_j += 2.*dof_rho;
						matrix->add(ind_i,ind_j, Q_drag*alpha_rho_g*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*alpha_rho_g*lumped);
//4.Zeile
						ind_i += dof_rho; ind_j = dof; ind_i_p = ind_i+dof_gas;
						matrix->add(ind_i,ind_j, Q_drag*((alpha_rho_v_x_p*alpha_rho_v_x_p+alpha_rho_v_y_p*alpha_rho_v_y_p)/alpha_rho_p)*lumped 
								-Q_tem *(alpha_rho_v_x_p/c_vg)*(-alpha_rho_e_g/(alpha_rho_g*alpha_rho_g) + (alpha_rho_v_x_g*alpha_rho_v_x_g+alpha_rho_v_y_g*alpha_rho_v_y_g)/(alpha_rho_g*alpha_rho_g*alpha_rho_g))*lumped);
						matrix->add(ind_i_p,ind_j, -Q_drag*((alpha_rho_v_x_p*alpha_rho_v_x_p+alpha_rho_v_y_p*alpha_rho_v_y_p)/alpha_rho_p)*lumped
						+ Q_tem *(alpha_rho_v_x_p/c_vg)*(-alpha_rho_e_g/(alpha_rho_g*alpha_rho_g) + (alpha_rho_v_x_g*alpha_rho_v_x_g+alpha_rho_v_y_g*alpha_rho_v_y_g)/(alpha_rho_g*alpha_rho_g*alpha_rho_g))*lumped);
 													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_v_x_p*lumped -Q_tem*(alpha_rho_v_x_p/c_vg)*(-alpha_rho_v_x_g)/(alpha_rho_g*alpha_rho_g)*lumped); 
						matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_v_x_p*lumped+ Q_tem*(alpha_rho_v_x_p/c_vg)*(-alpha_rho_v_x_g)/(alpha_rho_g*alpha_rho_g)*lumped); 
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*alpha_rho_v_y_p*lumped-Q_tem*(alpha_rho_v_x_p/c_vg)*(-alpha_rho_v_y_g)/(alpha_rho_g*alpha_rho_g)*lumped); 
						matrix->add(ind_i_p,ind_j, Q_drag*alpha_rho_v_y_p*lumped+ Q_tem*(alpha_rho_v_x_p/c_vg)*(-alpha_rho_v_y_g)/(alpha_rho_g*alpha_rho_g)*lumped); 
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_tem/alpha_rho_g*(alpha_rho_v_x_p/c_vg)*lumped);
						matrix->add(ind_i_p,ind_j, Q_tem/alpha_rho_g*(alpha_rho_v_x_p/c_vg)*lumped);

													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(alpha_rho_g*(alpha_rho_v_x_p*alpha_rho_v_x_p+alpha_rho_v_y_p*alpha_rho_v_y_p)/(alpha_rho_p*alpha_rho_p))*lumped
												 -Q_tem*(T_g-T_p+alpha_rho_e_p/alpha_rho_p-(v_x_p*v_x_p+v_y_p*v_y_p)/c_vp)*lumped);
						matrix->add(ind_i_p,ind_j, Q_drag*(alpha_rho_g*(alpha_rho_v_x_p*alpha_rho_v_x_p+alpha_rho_v_y_p*alpha_rho_v_y_p)/(alpha_rho_p*alpha_rho_p))*lumped
													+ Q_tem*(T_g-T_p+alpha_rho_e_p/alpha_rho_p-(v_x_p*v_x_p+v_y_p*v_y_p)/c_vp)*lumped);  
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(alpha_rho_v_x_g*alpha_rho_p-2.*alpha_rho_g*alpha_rho_v_x_p)/alpha_rho_p*lumped -Q_tem*(v_x_p/c_vp)*lumped);
						matrix->add(ind_i_p,ind_j, Q_drag*(alpha_rho_v_x_g*alpha_rho_p-2.*alpha_rho_g*alpha_rho_v_x_p)/alpha_rho_p*lumped+ Q_tem*(v_x_p/c_vp)*lumped);
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(alpha_rho_v_y_g*alpha_rho_p-2.*alpha_rho_g*alpha_rho_v_y_p)/alpha_rho_p*lumped -Q_tem*(v_y_p/c_vp)*lumped);  
						matrix->add(ind_i_p,ind_j, Q_drag*(alpha_rho_v_y_g*alpha_rho_p-2.*alpha_rho_g*alpha_rho_v_y_p)/alpha_rho_p*lumped+ Q_tem*(v_y_p/c_vp)*lumped);
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_tem*(-1./c_vp)*lumped); 
						matrix->add(ind_i_p,ind_j, Q_tem*(-1./c_vp)*lumped); 
					}
				
					}
				}
			 }
		}
		
	}
}


void multiply_matrices(CSCMatrix<double>* mat1,CSCMatrix<double>* mat2,CSCMatrix<double>* result)
{
printf("MatrixMulti"); 
int * Ap = mat2->get_Ap(); int *Ai = mat2->get_Ai(); double* Ax = mat2->get_Ax();

	result->zero(); int size = mat1->get_size(); double a, b;
        for(int i = 0; i < size; i++)        
			for(int k = 0; k<size; k++)
			{ double d = 0.;
          		for(int j = 0; j < Ap[i + 1] - Ap[i]; j++)
					if((a= mat1->get(k,Ai[Ap[i] + j]))!=0.)
						d+= Ax[Ap[i] + j]*a;
				
			if(std::fabs(d)>1e-8) result->add(k,i,d);
			}

       
        


	/*for(int i = 0; i<size; i++)	
		for(int k = 0; k<size; k++)
		{double d= 0.;
			for(int j = 0; j<size; j++)
			{	if((a= mat1->get(i,j))!=0.)				
					d += a*mat2->get(j,k);			
			}
			if(std::fabs(d)>1e-8) result->add(i,k,d);
		}*/
		

}






int main(int argc, char* argv[])
{


   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain_all.mesh", basemesh);
Element* e = NULL;Node* vn=NULL;


  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
	{ 
			basemesh->refine_all_elements();

		}

 	 mesh->copy(basemesh);


/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
   View::wait();*/


bool serendipity = true;
/*
		SpaceSharedPtr<double> space_rho_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_g(new H1Space<double>(mesh, P_INIT));


		SpaceSharedPtr<double> space_rho_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_p(new H1Space<double>(mesh, P_INIT));
*/

SpaceSharedPtr<double> space_rho_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e_g(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));

SpaceSharedPtr<double> space_rho_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_x_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_rho_v_y_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));	
SpaceSharedPtr<double> space_e_p(new L2_SEMI_CG_Space<double>(mesh, P_INIT, serendipity));


int dof_rho = space_rho_g->get_num_dofs();

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho_g, space_rho_v_x_g, space_rho_v_y_g, space_e_g,space_rho_p, space_rho_v_x_p, space_rho_v_y_p, space_e_p);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho_g(new CustomInitialCondition_rho(mesh, GAMMA));	
  MeshFunctionSharedPtr<double> init_rho_v_x_g(new  CustomInitialCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_v_y_g(new  ConstantSolution<double>(mesh,  0.));	
  MeshFunctionSharedPtr<double> init_rho_e_g(new CustomInitialCondition_e(mesh, GAMMA));

 MeshFunctionSharedPtr<double> init_rho_p(new CustomInitialCondition_rho(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_v_x_p(new  CustomInitialCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_v_y_p(new  ConstantSolution<double>(mesh,  0.));	
  MeshFunctionSharedPtr<double> init_rho_e_p(new CustomInitialCondition_e(mesh, GAMMA,true));


  MeshFunctionSharedPtr<double> bdry_rho_g(new BoundaryCondition_rho(mesh, GAMMA));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new  ConstantSolution<double>(mesh, 0.));	
  MeshFunctionSharedPtr<double> bdry_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));	

  MeshFunctionSharedPtr<double> bdry_rho_p(new BoundaryCondition_rho(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_p(new  BoundaryCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_p(new  ConstantSolution<double>(mesh, 0.));	
  MeshFunctionSharedPtr<double> bdry_rho_e_p(new BoundaryCondition_rho_e(mesh, GAMMA,true));	


  	MeshFunctionSharedPtr<double> prev_rho_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y_g(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_e_g(new Solution<double>);

  	MeshFunctionSharedPtr<double> prev_rho_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_x_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_v_y_p(new Solution<double>);
    MeshFunctionSharedPtr<double> prev_rho_e_p(new Solution<double>);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_g(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_g(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_p(prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_p(init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);


Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);


	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity

  ScalarView s1_g("rho", new WinGeom(0, 0, 400, 300));
 ScalarView s2_g("rho_v_x_g", new WinGeom(300, 0, 400, 300));
 ScalarView s3_g("rho_v_y_g", new WinGeom(300, 300, 400, 300));
ScalarView s4_g("rho_e_g", new WinGeom(300, 700, 400, 300));
  ScalarView pressure_view_g("Pressure-gas", new WinGeom(0, 300, 400, 300));
  ScalarView mach_view_g("mach-gas", new WinGeom(0, 750, 400, 300));


  ScalarView s1_p("rho_p", new WinGeom(700, 0, 400, 300));
 ScalarView s2_p("rho_v_x_p", new WinGeom(900, 0, 400, 300));
 ScalarView s3_p("rho_v_y_p", new WinGeom(900, 300, 400, 300));
ScalarView s4_p("rho_e_p", new WinGeom(900, 700, 400, 300));

  ScalarView pressure_view_p("Pressure_p", new WinGeom(700, 300, 400, 300));
  ScalarView mach_view_p("mach_p", new WinGeom(700, 750, 400, 300));
  ScalarView alpha_view_p("alpha_p", new WinGeom(700, 750, 400, 300));



s1_g.show(init_rho_g);
s1_p.show(init_rho_p);
s2_p.show(init_rho_v_x_p);
s3_p.show(init_rho_v_y_p);
s4_p.show(init_rho_e_p);


MeshFunctionSharedPtr<double> pressure_init_g(new PressureFilter(init_slns_g, GAMMA));
				pressure_view_g.show(pressure_init_g);
MeshFunctionSharedPtr<double> mach_init_g(new  MachNumberFilter(init_slns_g, GAMMA));
				mach_view_g.show(mach_init_g);





//View::wait(HERMES_WAIT_KEYPRESS);

//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(GAMMA);

NumericalFlux* num_flux =new LaxFriedrichsNumericalFlux(GAMMA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(GAMMA);



	EulerInterface wf_DG_init(GAMMA,mesh,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p,num_flux,euler_fluxes,riemann_invariants);
	EulerInterface wf_DG(GAMMA,mesh, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p,num_flux,euler_fluxes,riemann_invariants);


  DiscreteProblem<double> dp_DG_init(&wf_DG_init, spaces);
  DiscreteProblem<double> dp_DG(&wf_DG, spaces);



	EulerK wf_convection_init(GAMMA,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
	EulerK wf_convection(GAMMA, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
  EulerEquationsWeakForm_Mass wf_mass;

	EulerBoundary wf_bdry_init(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
	EulerBoundary wf_bdry(GAMMA, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

  EulerPenalty wf_penalty_init(SIGMA,density_particle,init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
EulerPenalty wf_penalty(SIGMA, density_particle, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

  EulerSource wf_source_init(density_particle,diameter,c_vg,c_vp,c_pg, Pr,mu, init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);
EulerSource wf_source(density_particle,diameter,c_vg,c_vp,c_pg, Pr,mu, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);


  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);

  DiscreteProblem<double> dp_conv_init(&wf_convection_init,spaces); 
  DiscreteProblem<double> dp_conv(&wf_convection, spaces); 

  DiscreteProblem<double> dp_bdry_init(&wf_bdry_init,spaces); 
  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 

  DiscreteProblem<double> dp_penalty_init(&wf_penalty_init,spaces); 
  DiscreteProblem<double> dp_penalty(&wf_penalty, spaces); 

  DiscreteProblem<double> dp_source_init(&wf_source_init,spaces); 
  DiscreteProblem<double> dp_source(&wf_source, spaces); 



  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;
		CSCMatrix<double> matrix2;    
	CSCMatrix<double>  mass_matrix;   
	CSCMatrix<double>  source_matrix;  
	CSCMatrix<double>  conv_matrix; 
	CSCMatrix<double>  penalty_matrix;
	CSCMatrix<double>  bdry_matrix;  
	CSCMatrix<double> dg_matrix; 

CSCMatrix<double> lumped_source_matrix;

CSCMatrix<double> proj_matrix;

    OGProjection<double> ogProjection;Lumped_Projection lumpedProjection;
		SimpleVector<double>  vec_dg(ndof);
		SimpleVector<double>  vec_rhs(ndof);
		SimpleVector<double>  vec_bdry(ndof);
		SimpleVector<double>  vec_conv(ndof);
		SimpleVector<double>  vec_source(ndof);
		SimpleVector<double>  vec_penalty(ndof);

		double* coeff_vec= new double[ndof];	
		double* coeff_vec_2 = new double[ndof];
		double* source_vec = new double[ndof];




//Projection of the initial condition
  //ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
	//		Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

lumpedProjection.project_lumped(spaces, init_slns, coeff_vec);
			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);



                      s1_g.show(prev_rho_g);                       
                        s2_g.show(prev_rho_v_x_g);
                        s3_g.show(prev_rho_v_y_g);
                        s4_g.show(prev_rho_e_g);
                        s1_p.show(prev_rho_p);
                        s2_p.show(prev_rho_v_x_p);
                        s3_p.show(prev_rho_v_y_p);
                        s4_p.show(prev_rho_e_p);


//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
double norm = 1000;
double norm_rel = 1000;
		Space<double>::assign_dofs(spaces);
		  dp_mass.assemble(&mass_matrix);
proj_matrix.create(mass_matrix.get_size(),mass_matrix.get_nnz(), mass_matrix.get_Ap(), mass_matrix.get_Ai(),mass_matrix.get_Ax());
mass_matrix.multiply_with_Scalar(1./time_step);
CSCMatrix<double> * lumped_matrix;
lumped_matrix = massLumping(&mass_matrix);

CSCMatrix<double> * lumped_matrix_proj;
lumped_matrix_proj = massLumping(&proj_matrix);
//Timestep loop
do
{	


 
	  
 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof); 		
 	// if(ts!=1){
			dp_conv.assemble(&conv_matrix, &vec_conv);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);
			dp_penalty.assemble(&penalty_matrix, &vec_penalty);
			dp_source.assemble(&source_matrix);
			//dp_source.assemble(&source_matrix, &vec_source);
				
		  //	dp_DG.assemble(&dg_matrix,&vec_dg);	
		/*}else{
			dp_conv_init.assemble(&conv_matrix, &vec_conv);
			dp_bdry_init.assemble(&bdry_matrix, &vec_bdry);
			dp_penalty_init.assemble(&penalty_matrix, &vec_penalty);
			dp_source_init.assemble(&source_matrix, &vec_source);
			
		 	//dp_DG_init.assemble(&dg_matrix,&vec_dg);	
		}*/

CSCMatrix<double>* diff = NULL;	
/*
			UMFPackLinearMatrixSolver<double> * proj_solver = new UMFPackLinearMatrixSolver<double> (lumped_matrix_proj,&vec_source);	
			try{
			 proj_solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}
		for(int i=0; i<ndof;i++)		
					source_vec[i] = proj_solver->get_sln_vector()[i];*/
assemble_vector_s(space_rho_g,prev_slns, ndof, dof_rho*4., dof_rho, coeff_vec_2, &source_matrix,lumped_matrix_proj);
lumped_matrix_proj->multiply_with_vector(coeff_vec_2, source_vec);





//matrix.create(dg_matrix.get_size(),dg_matrix.get_nnz(), dg_matrix.get_Ap(), dg_matrix.get_Ai(),dg_matrix.get_Ax());
//matrix.add_sparse_matrix(&source_matrix);
 diff =  artificialDiffusion(GAMMA,coeff_vec,spaces,&conv_matrix);

		matrix.create(source_matrix.get_size(),source_matrix.get_nnz(), source_matrix.get_Ap(), source_matrix.get_Ai(),source_matrix.get_Ax());
//matrix.create(conv_matrix.get_size(),conv_matrix.get_nnz(), conv_matrix.get_Ap(), conv_matrix.get_Ai(),conv_matrix.get_Ax());
//matrix.add_sparse_matrix(&dg_matrix);
		matrix.add_sparse_matrix(&bdry_matrix);
		matrix.add_sparse_matrix(&conv_matrix);
		matrix.add_sparse_matrix(diff);
		matrix.add_sparse_matrix(&penalty_matrix);
		matrix.multiply_with_Scalar(-theta); 
		matrix.add_sparse_matrix(lumped_matrix);  
		//matrix.add_sparse_matrix(&mass_matrix); 


matrix2.create(source_matrix.get_size(),source_matrix.get_nnz(), source_matrix.get_Ap(), source_matrix.get_Ai(),source_matrix.get_Ax());	
//matrix2.create(bdry_matrix.get_size(),bdry_matrix.get_nnz(), bdry_matrix.get_Ap(), bdry_matrix.get_Ai(),bdry_matrix.get_Ax());			
		matrix2.add_sparse_matrix(&bdry_matrix);
		matrix2.add_sparse_matrix(&penalty_matrix);
		matrix2.multiply_with_Scalar(-theta); 
		matrix2.add_sparse_matrix(lumped_matrix);  
		//matrix2.add_sparse_matrix(&mass_matrix); 
		
	//-------------rhs: ------------		
		matrix2.multiply_with_vector(coeff_vec, coeff_vec_2);
		//diff->multiply_with_vector(coeff_vec, coeff_vec_2);
		vec_rhs.zero(); 
		vec_rhs.add_vector(coeff_vec_2); 
		//vec_rhs.add_vector(&vec_dg); 
		vec_rhs.add_vector(&vec_bdry); 
		//vec_rhs.add_vector(&vec_conv); 
		vec_rhs.add_vector(&vec_penalty); 
		vec_rhs.add_vector(source_vec); 
//vec_rhs.add_vector(&vec_source); 


	//-------------------------solution of (M-theta(K+P+B+S)) (u(n+1)-u(n) = Sn +Ku(n) +Bn+Pn------------ 		
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (&matrix,&vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
	
		for(int i=0; i<ndof;i++)		
					coeff_vec[i] = solver->get_sln_vector()[i];							

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

	

			// Visualize the solution.
			Hermes::Mixins::Loggable::Static::info("Visualize"); 	
            MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns_g, GAMMA));
            sprintf(title, "Pressure gas: ts=%i",ts);
            pressure_view_g.set_title(title);
            pressure_g->reinit();
            pressure_view_g.show(pressure_g);

            MeshFunctionSharedPtr<double> mach_g(new MachNumberFilter(prev_slns_g, GAMMA));
            sprintf(title, "Mach gas: ts=%i",ts);
            mach_view_g.set_title(title);
            mach_g->reinit();
            mach_view_g.show(mach_g);


            MeshFunctionSharedPtr<double> alpha_p(new AlphaFilter(prev_slns_p, density_particle));
            sprintf(title, "alpha particle: ts=%i",ts);
            alpha_view_p.set_title(title);
            alpha_p->reinit();
            alpha_view_p.show(alpha_p);


            sprintf(title, "density__gas: ts=%i",ts);
            s1_g.set_title(title);
            s1_g.show(prev_rho_g);

			sprintf(title, "vel_x_g: ts=%i",ts);
            s2_g.set_title(title);
            s2_g.show(prev_rho_v_x_g);

			sprintf(title, "vel_y_g: ts=%i",ts);
            s3_g.set_title(title);
            s3_g.show(prev_rho_v_y_g);

			sprintf(title, "energy_g: ts=%i",ts);
            s4_g.set_title(title);
            s4_g.show(prev_rho_e_g);

			sprintf(title, "density_particle: ts=%i",ts);
            s1_p.set_title(title);
            s1_p.show(prev_rho_p);

			sprintf(title, "vel_x_particle: ts=%i",ts);
            s2_p.set_title(title);
            s2_p.show(prev_rho_v_x_p);

			sprintf(title, "vel_y_particle: ts=%i",ts);
            s3_p.set_title(title);
            s3_p.show(prev_rho_v_y_p);

			sprintf(title, "energy_particle: ts=%i",ts);
            s4_p.set_title(title);
            s4_p.show(prev_rho_e_p);


	//View::wait(HERMES_WAIT_KEYPRESS);




	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		matrix.free();
		if(diff!=NULL) delete diff;
 



}while (current_time < T_FINAL);




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  

/*
		MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, GAMMA));
		MeshFunctionSharedPtr<double> mach(new  MachNumberFilter(prev_slns, GAMMA));
 
Linearizer lin;	
				pressure->reinit();
				mach->reinit();
			  char filename[40];
			  sprintf(filename, "p-%i.vtk", ts );       
			lin.save_solution_vtk(pressure, filename, "pressure", true);
sprintf(filename, "m-%i.vtk", ts );      
			lin.save_solution_vtk(mach, filename, "mach", true);
sprintf(filename, "rho-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[0], filename, "density_particle", true);

*/




delete [] coeff_vec;
delete [] coeff_vec_2;


  // Wait for the view to be closed.
  View::wait();
  return 0;
}

