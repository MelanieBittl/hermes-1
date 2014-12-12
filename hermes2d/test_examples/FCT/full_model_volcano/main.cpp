#define HERMES_REPORT_ALL
#include "euler_util.h"
#include "definitions.h"
#include "initial_condition.h"
#include "euler_flux.h"
#include "lumped_projection.h"


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Solvers;

#include "mass_lumping.cpp"
#include "artificial_diffusion.cpp"
#include "fct.cpp"

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT =1;       						// Initial polynomial degree.
const double time_step = 1e-4;
const double T_FINAL = 2.;                       // Time interval length. 

const double theta = 1.;

// Equation parameters.  
 
     
// GAMMA.
const double GAMMA = 1.4; 

// Penalty Parameter.
double SIGMA = std::pow(10,3);

//Particle density_particle
const double density_particle = 2380.; 

const double diameter = 2e-4;//1e-5;		//10mu m
const double c_vg = 717.5;// Pelanti Diss: (gamma-1)=R/c_v
const double c_vp = 1300.;// Pelanti Diss
const double c_pg = 1004.5;// Pelanti Diss gamma = c_p/c_v
const double mu = 1e-5;// Pelanti Diss
const double kappa_g = 0.05;// Pelanti Diss
const double Pr =c_pg*mu/kappa_g;		// Pr= cp mu/kappa_g

const double g = 9.80665;// m/s2

 

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;






void assemble_source(SpaceSharedPtr<double> space,Hermes::vector<MeshFunctionSharedPtr<double> > slns, int dof_total, int dof_gas, int dof_rho, double* s, CSCMatrix<double>* matrix, CSCMatrix<double>* lumped_matrix)
{
	Element* e =NULL;
	AsmList<double> al;
	double x, y; int dof; double F_D_1, F_D_2, Q_T;
matrix->zero();
double kap = c_pg*mu/Pr;
bool visited[dof_gas]; 
for(int i= 0; i<dof_gas; i++) visited[i] = false; 

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
						
						double u_g_1 = slns[0]->get_pt_value(x,y,false,e)->val[0];
						double u_g_2 = slns[1]->get_pt_value(x,y,false,e)->val[0];
						double u_g_3 = slns[2]->get_pt_value(x,y,false,e)->val[0];
						double u_g_4 = slns[3]->get_pt_value(x,y,false,e)->val[0];
						
						double u_p_1 = slns[4]->get_pt_value(x,y,false,e)->val[0];
						double u_p_2 = slns[5]->get_pt_value(x,y,false,e)->val[0];
						double u_p_3 = slns[6]->get_pt_value(x,y,false,e)->val[0];
						double u_p_4 = slns[7]->get_pt_value(x,y,false,e)->val[0];
						
						

						double rho_g = u_g_1/alpha_g;  
						double rho_v_x_g = u_g_2/alpha_g; 
						double rho_v_y_g = u_g_3/alpha_g; 
						double rho_e_g = u_g_4/alpha_g;
						double v_x_g = rho_v_x_g/rho_g;
						double v_y_g = rho_v_y_g/rho_g;

						double rho_p = density_particle;  
						double rho_v_x_p = u_p_2/alpha_p; 
						double rho_v_y_p = u_p_3/alpha_p; 
						double rho_e_p = u_p_4/alpha_p;
						double v_x_p = rho_v_x_p/rho_p;
						double v_y_p = rho_v_y_p/rho_p;

						double v1_diff = (v_x_g - v_x_p);
						double v2_diff = (v_y_g - v_y_p);
						double v_diff_abs = std::sqrt(v1_diff*v1_diff+ v2_diff*v2_diff);

						double Re = rho_g*diameter*v_diff_abs/mu;						
						double C_D = 0.44;
						if(Re<1e-14) C_D =0;
						else if(Re<1000)
								C_D=24./Re*(1.+0.15*std::pow(Re,0.687));
						
						double Nu = 2.+0.65*std::sqrt(Re)*std::pow(Pr,1./3.);

						double T_g = 1./c_vg*(rho_e_g/rho_g-0.5*(v_x_g*v_x_g+v_y_g*v_y_g));
						double T_p = 1./c_vp*(rho_e_p/rho_p-0.5*(v_x_p*v_x_p+v_y_p*v_y_p));
						
						double Q_drag = 0.75*v_diff_abs*C_D/(diameter*alpha_g*density_particle);
						double Q_tem = Nu*6.*kap/(diameter*diameter*density_particle);

						F_D_1 = Q_drag* (u_g_2*u_p_1-u_p_2*u_g_1);//rho_g*alpha_p*C_D/diameter*v1_diff*v_diff_abs*0.75;
						F_D_2 = Q_drag*(u_g_3*u_p_1-u_p_3*u_g_1); //rho_g*alpha_p*C_D/diameter*v2_diff*v_diff_abs*0.75;
						Q_T = Q_tem*u_p_1*(T_g-T_p); //Nu*6.*kap/sqr(diameter)*alpha_p*(T_g-T_p);
						

						int ind = dof+dof_rho;
						s[ind] = -F_D_1;	s[ind+dof_gas] = F_D_1; ind += dof_rho;
						s[ind] = -F_D_2 + g*u_g_1;	s[ind+dof_gas] = F_D_2 + g*u_p_1; ind += dof_rho;
						s[ind] = -(F_D_1*v_x_p+F_D_2*v_y_p)-Q_T + g*u_g_3; s[ind+dof_gas] = (F_D_1*v_x_p+F_D_2*v_y_p)+Q_T+ g*u_p_3;	

//1. Zeile 0
//2.Zeile
						int ind_i = dof+dof_rho; int ind_j = dof; int ind_i_p = ind_i+dof_gas;
						double lumped = lumped_matrix->get(dof, dof);
						matrix->add(ind_i,ind_j, Q_drag*u_p_2*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*u_p_2*lumped); ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*u_p_1*lumped); matrix->add(ind_i_p,ind_j, Q_drag*u_p_1*lumped);  ind_j = dof+dof_gas;
						matrix->add(ind_i,ind_j, -Q_drag*u_g_2*lumped);matrix->add(ind_i_p,ind_j, Q_drag*u_g_2*lumped); ind_j += dof_rho;
						matrix->add(ind_i,ind_j, Q_drag*u_g_1*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*u_g_1*lumped); 
//3. Zeile
						ind_i += dof_rho; ind_j = dof; ind_i_p = ind_i+dof_gas;
						matrix->add(ind_i,ind_j, (Q_drag*u_p_3+g)*lumped);matrix->add(ind_i_p,ind_j, (-Q_drag*u_p_3)*lumped); ind_j += 2.*dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*u_p_1*lumped);matrix->add(ind_i_p,ind_j, Q_drag*u_p_1*lumped); ind_j = dof+dof_gas;
						matrix->add(ind_i,ind_j, -Q_drag*u_g_3*lumped);matrix->add(ind_i_p,ind_j, Q_drag*u_g_3*lumped+g*lumped); ind_j += 2.*dof_rho;
						matrix->add(ind_i,ind_j, Q_drag*u_g_1*lumped); matrix->add(ind_i_p,ind_j, -Q_drag*u_g_1*lumped);
//4.Zeile
						ind_i += dof_rho; ind_j = dof; ind_i_p = ind_i+dof_gas;
						matrix->add(ind_i,ind_j, Q_drag*(u_p_2*u_p_2+u_p_3*u_p_3)/u_p_1*lumped 
								-Q_tem *(u_p_1/c_vg)*(-u_g_4/(u_g_1*u_g_1)+ (u_g_2*u_g_2+u_g_3*u_g_3)/(u_g_1*u_g_1*u_g_1))*lumped);
						matrix->add(ind_i_p,ind_j, -Q_drag*(u_p_2*u_p_2+u_p_3*u_p_3)/u_p_1*lumped
						+ Q_tem *(u_p_1/c_vg)*(-u_g_4/(u_g_1*u_g_1)+ (u_g_2*u_g_2+u_g_3*u_g_3)/(u_g_1*u_g_1*u_g_1))*lumped);
						
 													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*u_p_2*lumped -Q_tem*(u_p_1/c_vg)*(-u_g_2)/(u_g_1*u_g_1)*lumped); 
						matrix->add(ind_i_p,ind_j, Q_drag*u_p_2*lumped+ Q_tem*(u_p_1/c_vg)*(-u_g_2)/(u_g_1*u_g_1)*lumped); 
						
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, (-Q_drag*u_p_3+g)*lumped-Q_tem*(u_p_1/c_vg)*(-u_g_3)/(u_g_1*u_g_1)*lumped); 
						matrix->add(ind_i_p,ind_j, (Q_drag*u_p_3)*lumped+ Q_tem*(u_p_1/c_vg)*(-u_g_3)/(u_g_1*u_g_1)*lumped); 
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_tem/u_g_1*(u_p_1/c_vg)*lumped);
						matrix->add(ind_i_p,ind_j, Q_tem/u_g_1*(u_p_1/c_vg)*lumped);

													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(u_g_1*(u_p_2*u_p_2+u_p_3*u_p_3)/(u_p_1*u_p_1))*lumped
												 -Q_tem*(T_g-T_p+u_p_4/(c_vp*u_p_1)-(u_p_2*u_p_2+u_p_3*u_p_3)/(c_vp*u_p_1*u_p_1))*lumped);
						matrix->add(ind_i_p,ind_j, Q_drag*(u_g_1*(u_p_2*u_p_2+u_p_3*u_p_3)/(u_p_1*u_p_1))*lumped
													+ Q_tem*(T_g-T_p+u_p_4/(c_vp*u_p_1)-(u_p_2*u_p_2+u_p_3*u_p_3)/(c_vp*u_p_1*u_p_1))*lumped);  
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(u_g_2*u_p_1-2.*u_g_1*u_p_2)/u_p_1*lumped -Q_tem*u_p_2/(c_vp*u_p_1)*lumped);
						matrix->add(ind_i_p,ind_j, Q_drag*(u_g_2*u_p_1-2.*u_g_1*u_p_2)/u_p_1*lumped+ Q_tem*u_p_2/(c_vp*u_p_1)*lumped);
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, -Q_drag*(u_g_3*u_p_1-2.*u_g_1*u_p_3)/u_p_1*lumped -Q_tem*u_p_3/(c_vp*u_p_1)*lumped);  
						matrix->add(ind_i_p,ind_j, Q_drag*(u_g_3*u_p_1-2.*u_g_1*u_p_3)/u_p_1*lumped+ Q_tem*u_p_3/(c_vp*u_p_1)*lumped+g*lumped);
													 ind_j += dof_rho;
						matrix->add(ind_i,ind_j, Q_tem/c_vp*lumped); 
						matrix->add(ind_i_p,ind_j, -Q_tem/c_vp*lumped); 
					}
				
					}
				}
			 }
		}
		
	}
}


void assemble_vector_s(SpaceSharedPtr<double> space,Hermes::vector<MeshFunctionSharedPtr<double> > slns, int dof_total, int dof_gas, int dof_rho, double* s)
{
	Element* e =NULL;
	AsmList<double> al;
	double x, y; int dof; double F_D_1, F_D_2, Q_T;

double kap = c_pg*mu/Pr;
bool visited[dof_gas]; 
for(int i= 0; i<dof_gas; i++) visited[i] = false; 

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
						
						double u_g_1 = slns[0]->get_pt_value(x,y,false,e)->val[0];
						double u_g_2 = slns[1]->get_pt_value(x,y,false,e)->val[0];
						double u_g_3 = slns[2]->get_pt_value(x,y,false,e)->val[0];
						double u_g_4 = slns[3]->get_pt_value(x,y,false,e)->val[0];
						
						double u_p_1 = slns[4]->get_pt_value(x,y,false,e)->val[0];
						double u_p_2 = slns[5]->get_pt_value(x,y,false,e)->val[0];
						double u_p_3 = slns[6]->get_pt_value(x,y,false,e)->val[0];
						double u_p_4 = slns[7]->get_pt_value(x,y,false,e)->val[0];				
						

						double rho_g = u_g_1/alpha_g;  
						double rho_v_x_g = u_g_2/alpha_g; 
						double rho_v_y_g = u_g_3/alpha_g; 
						double rho_e_g = u_g_4/alpha_g;
						double v_x_g = rho_v_x_g/rho_g;
						double v_y_g = rho_v_y_g/rho_g;

						double rho_p = density_particle;  
						double rho_v_x_p = u_p_2/alpha_p; 
						double rho_v_y_p = u_p_3/alpha_p; 
						double rho_e_p = u_p_4/alpha_p;
						double v_x_p = rho_v_x_p/rho_p;
						double v_y_p = rho_v_y_p/rho_p;

						double v1_diff = (v_x_g - v_x_p);
						double v2_diff = (v_y_g - v_y_p);
						double v_diff_abs = std::sqrt(v1_diff*v1_diff+ v2_diff*v2_diff);

						double Re = rho_g*diameter*v_diff_abs/mu;						
						double C_D = 0.44;
						if(Re<1e-14) C_D =0;
						else if(Re<1000)
								C_D=24./Re*(1.+0.15*std::pow(Re,0.687));
						
						double Nu = 2.+0.65*std::sqrt(Re)*std::pow(Pr,1./3.);

						double T_g = 1./c_vg*(rho_e_g/rho_g-0.5*(v_x_g*v_x_g+v_y_g*v_y_g));
						double T_p = 1./c_vp*(rho_e_p/rho_p-0.5*(v_x_p*v_x_p+v_y_p*v_y_p));
						
						double Q_drag = 0.75*v_diff_abs*C_D/(diameter*alpha_g*density_particle);
						double Q_tem = Nu*6.*kap/(diameter*diameter*density_particle);

						F_D_1 = Q_drag* (u_g_2*u_p_1-u_p_2*u_g_1);
						F_D_2 = Q_drag*(u_g_3*u_p_1-u_p_3*u_g_1); 
						Q_T = Q_tem*u_p_1*(T_g-T_p); 

					
						int ind = dof+dof_rho;
						s[ind] = -F_D_1;	s[ind+dof_gas] = F_D_1; ind += dof_rho;
						s[ind] = -F_D_2 + g*u_g_1;	s[ind+dof_gas] = F_D_2 + g*u_p_1; ind += dof_rho;
						s[ind] = -(F_D_1*v_x_p+F_D_2*v_y_p)-Q_T + g*u_g_3; s[ind+dof_gas] = (F_D_1*v_x_p+F_D_2*v_y_p)+Q_T+ g*u_p_3;	

					}
				
					}
				}
			 }
		}
		
	}
}





int main(int argc, char* argv[])
{

   // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain3.mesh", basemesh);
Element* e = NULL;Node* vn=NULL;

 //domain3
basemesh->refine_in_area("C",1,2);		
  // Perform initial mesh refinements (optional).
 for (int i=0; i < INIT_REF_NUM; i++)
	{ 			//basemesh->refine_all_elements(1);
			//basemesh->refine_all_elements();
basemesh->refine_in_area("B");
basemesh->refine_in_area("C");
basemesh->refine_in_area("D",1,1);
basemesh->refine_in_area("E");

		}
basemesh->refine_in_area("D",2);
basemesh->refine_in_area("E");
//basemesh->refine_in_area("A",2,2);
basemesh->refine_in_area("B",2,1);	
basemesh->refine_in_area("C",2,1);	

/*
//domain
 for (int i=0; i < (INIT_REF_NUM+1); i++)
	{ 	
basemesh->refine_in_area("B");
basemesh->refine_in_area("D",1,1);

		}*/
//basemesh->refine_towards_boundary("in");



 	 mesh->copy(basemesh);




/*
   MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
 meshview.show(mesh);
  
   View::wait();

*/


	SpaceSharedPtr<double> space_rho_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_g(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_g(new H1Space<double>(mesh, P_INIT));

		SpaceSharedPtr<double> space_rho_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_x_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_rho_v_y_p(new H1Space<double>(mesh, P_INIT));	
		SpaceSharedPtr<double> space_e_p(new H1Space<double>(mesh, P_INIT));

int dof_rho = space_rho_g->get_num_dofs();

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho_g, space_rho_v_x_g, space_rho_v_y_g, space_e_g,space_rho_p, space_rho_v_x_p, space_rho_v_y_p, space_e_p);
	Space<double>::assign_dofs(spaces);
   int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d \n", ndof);

  // Initialize solutions, set initial conditions.
 MeshFunctionSharedPtr<double> init_rho_g(new CustomInitialCondition_rho(mesh, GAMMA));	
  MeshFunctionSharedPtr<double> init_rho_v_x_g(new  CustomInitialCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_v_y_g(new  CustomInitialCondition_rho_v_y(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> init_rho_e_g(new CustomInitialCondition_e(mesh, GAMMA));

 MeshFunctionSharedPtr<double> init_rho_p(new CustomInitialCondition_rho(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_v_x_p(new  CustomInitialCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_v_y_p(new  CustomInitialCondition_rho_v_y(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> init_rho_e_p(new CustomInitialCondition_e(mesh, GAMMA,true));

  

  MeshFunctionSharedPtr<double> bdry_rho_g(new BoundaryCondition_rho(mesh, GAMMA));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_g(new  BoundaryCondition_rho_v_x(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_g(new   BoundaryCondition_rho_v_y(mesh, GAMMA) );	
  MeshFunctionSharedPtr<double> bdry_rho_e_g(new BoundaryCondition_rho_e(mesh, GAMMA));	

  MeshFunctionSharedPtr<double> bdry_rho_p(new BoundaryCondition_rho(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_v_x_p(new  BoundaryCondition_rho_v_x(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_v_y_p(new  BoundaryCondition_rho_v_y(mesh, GAMMA,true));	
  MeshFunctionSharedPtr<double> bdry_rho_e_p(new BoundaryCondition_rho_e(mesh, GAMMA,true));	
 
 
   	MeshFunctionSharedPtr<double> low_rho_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_x_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_v_y_g(new Solution<double>);
    MeshFunctionSharedPtr<double> low_rho_e_g(new Solution<double>);

  	MeshFunctionSharedPtr<double> low_rho_p(new Solution<double>);
	MeshFunctionSharedPtr<double> low_rho_v_x_p(new Solution<double>);
	MeshFunctionSharedPtr<double> low_rho_v_y_p(new Solution<double>);
	MeshFunctionSharedPtr<double> low_rho_e_p(new Solution<double>);

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
	Hermes::vector<MeshFunctionSharedPtr<double> > low_slns_g(low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g);

	Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns_p(prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
	Hermes::vector<MeshFunctionSharedPtr<double> > init_slns_p(init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);		
				Hermes::vector<MeshFunctionSharedPtr<double> > low_slns_p(low_rho_p, low_rho_v_x_p, low_rho_v_y_p, low_rho_e_p);


Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);

Hermes::vector<MeshFunctionSharedPtr<double> > init_slns(init_rho_g, init_rho_v_x_g, init_rho_v_y_g, init_rho_e_g,init_rho_p, init_rho_v_x_p, init_rho_v_y_p, init_rho_e_p);

Hermes::vector<MeshFunctionSharedPtr<double> > low_slns(low_rho_g, low_rho_v_x_g, low_rho_v_y_g, low_rho_e_g,low_rho_p,low_rho_v_x_p,low_rho_v_y_p,low_rho_e_p);


	Hermes::vector<NormType> norms_l2(HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM,HERMES_L2_NORM);


 //--------- Visualization of pressure & velocity

  ScalarView s1_g("rho", new WinGeom(0, 0, 400, 300));
 ScalarView s2_g("rho_v_x_g", new WinGeom(300, 0, 400, 300));
 ScalarView s3_g("rho_v_y_g", new WinGeom(300, 300, 400, 300));
ScalarView s4_g("rho_e_g", new WinGeom(300, 600, 400, 300));
  ScalarView pressure_view_g("Pressure-gas", new WinGeom(0, 300, 400, 300));
  ScalarView mach_view_g("mach-gas", new WinGeom(0, 600, 400, 300));
  ScalarView temp_view_g("Temperature_g", new WinGeom(0, 900, 400, 300));

  ScalarView s1_p("rho_p", new WinGeom(700, 0, 400, 300));
 ScalarView s2_p("rho_v_x_p", new WinGeom(900, 0, 400, 300));
 ScalarView s3_p("rho_v_y_p", new WinGeom(900, 300, 400, 300));
ScalarView s4_p("rho_e_p", new WinGeom(900, 600, 400, 300));

  ScalarView temp_view_p("Temperature_p", new WinGeom(700, 300, 400, 300));
  ScalarView mach_view_p("mach_p", new WinGeom(700, 750, 400, 300));
  ScalarView alpha_view_p("alpha_p", new WinGeom(700, 0, 400, 300));


/*
s1_g.show(init_rho_g);
                  s2_g.show(init_rho_v_x_g);
                        s3_g.show(init_rho_v_y_g);
                        s4_g.show(init_rho_e_g);
s1_p.show(init_rho_p);
s2_p.show(init_rho_v_x_p);
s3_p.show(init_rho_v_y_p);
s4_p.show(init_rho_e_p);
*/

MeshFunctionSharedPtr<double> pressure_init_g(new PressureFilter(init_slns_g, GAMMA));
			//	pressure_view_g.show(pressure_init_g);
MeshFunctionSharedPtr<double> mach_init_g(new  MachNumberFilter(init_slns_g, GAMMA));
			//	mach_view_g.show(mach_init_g);





//View::wait(HERMES_WAIT_KEYPRESS);

//------------

	EulerFluxes* euler_fluxes = new EulerFluxes(GAMMA);

	RiemannInvariants* riemann_invariants = new RiemannInvariants(GAMMA);

	EulerK wf_convection(GAMMA, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);
  EulerEquationsWeakForm_Mass wf_mass;


	EulerBoundary wf_bdry(GAMMA,SIGMA,density_particle, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);	
	
	EulerBoundary wf_bdry_low(GAMMA,SIGMA,density_particle, bdry_rho_g,bdry_rho_v_x_g,bdry_rho_v_y_g, bdry_rho_e_g, bdry_rho_p,bdry_rho_v_x_p,bdry_rho_v_y_p, bdry_rho_e_p, low_rho_g, low_rho_v_x_g,low_rho_v_y_g, low_rho_e_g,low_rho_p, low_rho_v_x_p, low_rho_v_y_p, low_rho_e_p);

EulerSource wf_source(density_particle,diameter,c_vg,c_vp,c_pg, Pr,mu, prev_rho_g, prev_rho_v_x_g, prev_rho_v_y_g, prev_rho_e_g,prev_rho_p, prev_rho_v_x_p, prev_rho_v_y_p, prev_rho_e_p);



  // Initialize the FE problem.
  DiscreteProblem<double> dp_mass(&wf_mass,spaces);
 
  DiscreteProblem<double> dp_conv(&wf_convection, spaces); 

  DiscreteProblem<double> dp_bdry(&wf_bdry, spaces); 
  DiscreteProblem<double> dp_bdry_low(&wf_bdry_low, spaces); 

  DiscreteProblem<double> dp_source(&wf_source, spaces);   





  // Set up the solver, matrix, and rhs according to the solver selection. 
	CSCMatrix<double> matrix;
		CSCMatrix<double> matrix2;    
	CSCMatrix<double>  mass_matrix;   
	CSCMatrix<double>  source_matrix;  

	CSCMatrix<double>  bdry_matrix; 
	CSCMatrix<double> matrixL; 



CSCMatrix<double> proj_matrix;

    OGProjection<double> ogProjection;Lumped_Projection lumpedProjection;
		SimpleVector<double>  vec_rhs(ndof);
		SimpleVector<double>  vec_bdry(ndof);
		SimpleVector<double>  vec_conv(ndof);
		SimpleVector<double>  vec_source(ndof);

	SimpleVector<double> * vec_res = new SimpleVector<double> (ndof);

		double* coeff_vec= new double[ndof];	
		double* coeff_vec_2 = new double[ndof];
		double* source_vec = new double[ndof];
		
		double* u_L= new double[ndof];
		double* flux_vec = new double[ndof];
		double* P_plus = new double[ndof]; double* P_minus = new double[ndof];
		double* Q_plus = new double[ndof]; double* Q_minus = new double[ndof];
		double* R_plus = new double[ndof]; double* R_minus = new double[ndof];	




//Projection of the initial condition
  //ogProjection.project_global(spaces,init_slns, coeff_vec, norms_l2 );
		//Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);

lumpedProjection.project_lumped(spaces, init_slns, coeff_vec);
		Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);
		
		
		/*	MeshFunctionSharedPtr<double> pressure_g(new PressureFilter(prev_slns_g, GAMMA));
			MeshFunctionSharedPtr<double> mach_g(new MachNumberFilter(prev_slns_g, GAMMA));
			MeshFunctionSharedPtr<double> alpha_p(new AlphaFilter(prev_slns_p, density_particle));
			MeshFunctionSharedPtr<double> temp_p(new TempFilter(prev_slns_p, c_vp));
			MeshFunctionSharedPtr<double> temp_g(new TempFilter(prev_slns_g, c_vg));*/

Linearizer lin;	

               /*     s1_g.show(prev_rho_g);                       
                        s2_g.show(prev_rho_v_x_g);
                        s3_g.show(prev_rho_v_y_g);
                        s4_g.show(prev_rho_e_g);
                        s1_p.show(prev_rho_p);
                        s2_p.show(prev_rho_v_x_p);
                        s3_p.show(prev_rho_v_y_p);
                        s4_p.show(prev_rho_e_p);*/


//View::wait(HERMES_WAIT_KEYPRESS);

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];	
	  char filename[40];
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

double residual = 1000000;
dp_source.assemble(&matrix); 
//Timestep loop
do
{	

 Hermes::Mixins::Loggable::Static::info("Time step %d,  time %3.5f, ndofs=%i", ts, current_time, ndof); 		
			dp_conv.assemble(&matrixL);
			dp_bdry.assemble(&bdry_matrix, &vec_bdry);

CSCMatrix<double>* diff = NULL;	
assemble_source(space_rho_g,prev_slns, ndof, dof_rho*4., dof_rho, coeff_vec_2, &matrix,lumped_matrix_proj);
lumped_matrix_proj->multiply_with_vector(coeff_vec_2, source_vec);

diff =  artificialDiffusion(GAMMA,coeff_vec,spaces,&matrixL);
matrixL.add_sparse_matrix(diff);
matrix.add_sparse_matrix(&bdry_matrix); // source + bdry

	//matrix.create(matrix2.get_size(),matrix2.get_nnz(), matrix2.get_Ap(),matrix2.get_Ai(),matrix2.get_Ax());
		matrix.add_sparse_matrix(&matrixL);
		matrix.multiply_with_Scalar(-theta); 
		matrix.add_sparse_matrix(lumped_matrix); 

		//matrix2.multiply_with_Scalar(-theta); 
		//matrix2.add_sparse_matrix(lumped_matrix);  


	//-------------rhs: ------------		
		matrixL.multiply_with_vector(coeff_vec, coeff_vec_2);//KU+DU
		vec_rhs.zero(); 
		vec_rhs.add_vector(coeff_vec_2); 
		vec_rhs.add_vector(&vec_bdry); 
		vec_rhs.add_vector(source_vec); 





	//-------------------------solution of (M-theta(K+P+B+S)) (u(n+1)-u(n) = Sn +Ku(n) +Bn+Pn------------ 		
			UMFPackLinearMatrixSolver<double> * solver = new UMFPackLinearMatrixSolver<double> (&matrix,&vec_rhs);	
			try{
			 solver->solve();
			}catch(Hermes::Exceptions::Exception e){
				e.print_msg();
			}	
	
		for(int i=0; i<ndof;i++)	
				coeff_vec[i]+= solver->get_sln_vector()[i];		
	
		
/*
	for(int i=0; i<ndof;i++)	
								u_L[i] = solver->get_sln_vector()[i] + coeff_vec[i];	
	
	Solution<double>::vector_to_solutions(u_L, spaces, low_slns);
		dp_bdry_low.assemble(&vec_bdry);		
			assemble_vector_s(space_rho_g,low_slns, ndof, dof_rho*4., dof_rho, coeff_vec_2);
lumped_matrix_proj->multiply_with_vector(coeff_vec_2, source_vec);


antidiffusiveFlux(&mass_matrix,lumped_matrix,diff,&matrixL, &vec_bdry,source_vec, u_L, flux_vec, P_plus,P_minus, Q_plus, Q_minus,  R_plus, R_minus,dof_rho,time_step,GAMMA );
	
		for(int i=0; i<ndof;i++)		
					coeff_vec[i] = u_L[i] + flux_vec[i]*time_step/lumped_matrix->get(i,i);
*/

			Solution<double>::vector_to_solutions(coeff_vec, spaces, prev_slns);				
			

	
	// Visualize the solution.
/*Hermes::Mixins::Loggable::Static::info("Visualize"); 									
					
					
							sprintf(title, "density_particle: ts=%i",ts);
            s1_p.set_title(title);
            s1_p.show(prev_rho_p);
					
            sprintf(title, "Pressure gas: ts=%i",ts);
            pressure_view_g.set_title(title);
            pressure_g->reinit();
            pressure_view_g.show(pressure_g);

            sprintf(title, "Mach gas: ts=%i",ts);
            mach_view_g.set_title(title);
            mach_g->reinit();
            mach_view_g.show(mach_g);

            sprintf(title, "alpha particle: ts=%i",ts);
            alpha_view_p.set_title(title);
            alpha_p->reinit();
            alpha_view_p.show(alpha_p);
				

            sprintf(title, "Temp gas: ts=%i",ts);
            temp_view_g.set_title(title);
            temp_g->reinit();
            temp_view_g.show(temp_g);

          /*  sprintf(title, "Temp particle: ts=%i",ts);
            temp_view_p.set_title(title);
            temp_p->reinit();
            temp_view_p.show(temp_p);

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
            s4_p.show(prev_rho_e_p);*/


	//View::wait(HERMES_WAIT_KEYPRESS);

	if(ts%200 ==0)
{

			/*	alpha_p->reinit();
sprintf(filename, "alpha_p-%i.vtk", ts );  
			lin.save_solution_vtk(alpha_p, filename, "alpha_particle", false);*/
			sprintf(filename, "rho_p-%i.vtk", ts );  
			lin.save_solution_vtk(prev_rho_p, filename, "alpha_particle", false);

}

/*	if(ts%1000 ==0)
{
				pressure_g->reinit();
				mach_g->reinit();
				      temp_g->reinit();
			
			  sprintf(filename, "p_gas_%i.vtk", ts );       
			lin.save_solution_vtk(pressure_g, filename, "pressure", false);
						  sprintf(filename, "t_gas_%i.vtk", ts );       
			lin.save_solution_vtk(temp_g, filename, "temperature_gas", false);
sprintf(filename, "m_gas_%i.vtk", ts );      
			lin.save_solution_vtk(mach_g, filename, "mach", false);
sprintf(filename, "rho_g-%i.vtk", ts );  
			lin.save_solution_vtk(prev_slns[0], filename, "density_gas", false);


}*/



	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;

		delete solver;
		//matrix.free();
		if(diff!=NULL) delete diff;
 



}while (current_time < T_FINAL);




Hermes::Mixins::Loggable::Static::info("end_time %3.5f",current_time); 	  
/*
				pressure_g->reinit();
				mach_g->reinit();
				alpha_p->reinit();
				temp_g->reinit();
			
    
			lin.save_solution_vtk(pressure_g, "p_gas_end.vtk", "pressure", false);    
			lin.save_solution_vtk(mach_g, "m_gas_end.vtk", "mach", false); 
			lin.save_solution_vtk(prev_slns[0], "rho_gas_end.vtk", "density_gas", false);
			lin.save_solution_vtk(alpha_p, "alpha_p_end.vtk", "alpha_particle", false);
			lin.save_solution_vtk(temp_g, "temp_gas_end.vtk", "temperature_gas", false);
*/

delete [] coeff_vec;
delete [] coeff_vec_2;


delete [] source_vec;
		
delete [] u_L;
delete [] flux_vec;
delete [] P_plus; delete [] P_minus ;
delete [] Q_plus; delete [] Q_minus ;
delete [] R_plus; delete [] R_minus ;	
delete vec_res;



  // Wait for the view to be closed.
  View::wait();
  return 0;
}

