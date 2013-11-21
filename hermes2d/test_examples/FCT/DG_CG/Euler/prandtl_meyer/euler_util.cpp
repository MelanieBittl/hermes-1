#include "euler_util.h"
#include "limits.h"
#include <limits>

// Calculates energy from other quantities.
double QuantityCalculator::calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa)
{
  double to_return = pressure/(kappa - 1.0) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / (2.0*rho);
 // if(std::abs(to_return) < 1E-12 || to_return < 0.0)
   // return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

// Calculates pressure from other quantities.
double QuantityCalculator::calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = (kappa - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2.0*rho));
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

// Calculates speed of sound.
double QuantityCalculator::calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa)
{
  double to_return = std::sqrt(kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy, kappa) / rho);
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
  if(to_return < 0.0)
    return 0.0;
  return to_return;
}

  // Calculates enthalpy.
double QuantityCalculator::enthalpy(double rho, double rho_v_x, double rho_v_y, double energy, double kappa){

double to_return= (energy+ QuantityCalculator::calc_pressure(rho, rho_v_x, rho_v_y,energy,kappa))/rho;
  //if(std::abs(to_return) < 1E-12 || to_return < 0.0)
    //return 1E-12;
 if(to_return < 0.0)
   return 0.0;

return to_return;
}

//----------------------------------------------------------------------------------
//----------------------Filters.--------------------------------------
//----------------------------------------------------------------------------------
void MachNumberFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt((values.at(1)[i] / values.at(0)[i])*(values.at(1)[i] / values.at(0)[i]) + (values.at(2)[i] / values.at(0)[i])*(values.at(2)[i] / values.at(0)[i]))
    / std::sqrt(kappa * QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / values.at(0)[i]);
}


void PressureFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = (kappa - 1.) * (values.at(3)[i] - (values.at(1)[i]*values.at(1)[i] + values.at(2)[i]*values.at(2)[i])/(2*values.at(0)[i]));
}

void VelocityFilter_x::filter_fn(int n, Hermes::vector<double*> values, double* result)
{ 
  for (int i = 0; i < n; i++)
    result[i] = values.at(1)[i]/(values.at(0)[i]);
}
void VelocityFilter_y::filter_fn(int n, Hermes::vector<double*> values, double* result)
{ 
  for (int i = 0; i < n; i++)
    result[i] = values.at(2)[i]/(values.at(0)[i]);
}

void RadiusVelocityFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
{
  for (int i = 0; i < n; i++)
    result[i] = std::sqrt( values.at(1)[i]*values.at(1)[i]+values.at(2)[i]*values.at(2)[i])/values.at(0)[i];
}

void EntropyFilter::filter_fn(int n, Hermes::vector<double*> values, double* result) 
{
  for (int i = 0; i < n; i++)
    for (int i = 0; i < n; i++)
      result[i] = std::log((QuantityCalculator::calc_pressure(values.at(0)[i], values.at(1)[i], values.at(2)[i], values.at(3)[i], kappa) / p_ext)
      / Hermes::pow((values.at(0)[i] / rho_ext), kappa));
}


//-----------------------------------------------------
//-------------------RiemannInvariants-----------------------
//-----------------------------------------------------
double RiemannInvariants::get_w1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - 2*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa)/(kappa-1));
}

double RiemannInvariants::get_w2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return QuantityCalculator::calc_pressure(rho, rho_v_x,rho_v_y, rho_energy, kappa)/Hermes::pow(rho,kappa);
}

double RiemannInvariants::get_w3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double t_x, double t_y){
		return  (rho_v_x*t_x/rho + rho_v_y*t_y/rho);
}

double RiemannInvariants::get_w4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +2*QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa)/(kappa-1));
}

bool RiemannInvariants::solid_wall(double rho, double rho_v_x, double rho_v_y, double n_x, double n_y){
	if(fabs(rho_v_x*n_x/rho + rho_v_y*n_y/rho)<1E-12) 
				return true;
	else false;

}


double RiemannInvariants::get_ev1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho - QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa));
}

double RiemannInvariants::get_ev2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
		return (rho_v_x*n_x/rho + rho_v_y*n_y/rho);
}

double RiemannInvariants::get_ev4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y){
	return (rho_v_x*n_x/rho + rho_v_y*n_y/rho +QuantityCalculator::calc_sound_speed(rho, rho_v_x,rho_v_y, rho_energy, kappa));
}

double RiemannInvariants::get_pressure(double w_1, double w_2, double w_3, double w_4){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return (c*c*rho/kappa);
}

double RiemannInvariants::get_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x){
	return (w_1+w_4)*n_x/2.+w_3*t_x;
}
double RiemannInvariants::get_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y){
	return (w_1+w_4)*n_y/2.+w_3*t_y;

}
double RiemannInvariants::get_speed_sound(double w_1, double w_4){
									return (kappa-1.)*(w_4-w_1)/4.;
}

double RiemannInvariants::get_rho(double w_1, double w_2, double w_3, double w_4){
	double c = RiemannInvariants::get_speed_sound(w_1,  w_4);
	return  Hermes::pow(c*c/(kappa*w_2),(1/(kappa-1)));
}
double RiemannInvariants::get_rho_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double v_x = RiemannInvariants::get_v_x(w_1,w_2, w_3, w_4, n_x, t_x);
	return v_x*rho;

}
double RiemannInvariants::get_rho_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double v_y = RiemannInvariants::get_v_y(w_1,w_2, w_3, w_4, n_y, t_y);
		return v_y*rho;

}
double RiemannInvariants::get_energy(double w_1, double w_2, double w_3, double w_4, double n_x, double n_y, double t_x, double t_y){
	double rho = RiemannInvariants::get_rho(w_1,w_2, w_3, w_4);
	double p = RiemannInvariants::get_pressure(w_1,w_2, w_3, w_4);
	double v_x = RiemannInvariants::get_v_x(w_1,w_2, w_3, w_4, n_x, t_x);
	double v_y = RiemannInvariants::get_v_y(w_1,w_2, w_3, w_4, n_y, t_y);
	double v_norm = v_x*v_x+v_y*v_y;	
	return (p/(kappa-1.)+rho*v_norm/2.);

}


void RiemannInvariants::get_ghost_state(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid)
{

		double w_1,w_2,w_3,w_4;	

	if( solid==true){
		ghost_state[0] = rho;
		ghost_state[3] = rho_energy;
		ghost_state[1] = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
		ghost_state[2] = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
	}else if(bdry==1){ //supersonic outlet{
		ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=		rho_v_y;
		ghost_state[3]= rho_energy;
	}else if(bdry==2){//supersonic inlet
		ghost_state[0] =  rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=		rho_v_y_ext;
		ghost_state[3]= rho_energy_ext;
	
	}
}

int RiemannInvariants::get_bdry_info(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid){

	double lambda_1 = RiemannInvariants::get_ev1(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_2_3 = RiemannInvariants::get_ev2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_4 = RiemannInvariants::get_ev4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
		double w_1,w_2,w_3,w_4;	

	if( solid==true){
		ghost_state[0] = rho;
		ghost_state[3] = rho_energy;
		ghost_state[1] = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
		ghost_state[2] = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
			return 0;
	}else if((lambda_1>=0)&&(lambda_2_3>=0)&&(lambda_4>=0)){ //supersonic outlet{
		ghost_state[0] =  rho;
		ghost_state[1] = rho_v_x;
		ghost_state[2]=		rho_v_y;
		ghost_state[3]= rho_energy;
			return 1;
	}else if((lambda_1<0)&&(lambda_2_3<0)&&(lambda_4<0)){//supersonic inlet
		ghost_state[0] =  rho_ext;
		ghost_state[1] = rho_v_x_ext;
		ghost_state[2]=		rho_v_y_ext;
		ghost_state[3]= rho_energy_ext;
			return 2;			
	}else return 5;
}

void RiemannInvariants::get_du_du(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext,int bdry, int entry_j, double* dudu){


	double c = QuantityCalculator::calc_sound_speed(rho, rho_v_x, rho_v_y, rho_energy,kappa);
	double c_ext = QuantityCalculator::calc_sound_speed(rho_ext, rho_v_x_ext, rho_v_y_ext, rho_energy_ext,kappa);

 if(bdry ==0){ //solid wall
				if(entry_j==0){
						dudu[0] = 1.;
						dudu[1] = 0.;
						dudu[2] = 0.;
						dudu[3] = 0.;
				}else if(entry_j==1){
						dudu[0] = 0.;
						dudu[1] = (n_y*n_y-n_x*n_x);
						dudu[2] = (-2*n_x*n_y);
						dudu[3] = 0.;
				}else if(entry_j==2){
						dudu[0] = 0.;
						dudu[1] = (-2*n_x*n_y);
						dudu[2] = (n_x*n_x-n_y*n_y);
						dudu[3] = 0.;
				}else if(entry_j==3){
						dudu[0] = 0.;
						dudu[1] = 0.;
						dudu[2] = 0.;
						dudu[3] = 1.;
				}else  throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no entry_j!");


	}else   throw Hermes::Exceptions::Exception("RiemannInvariants::get_du_du: no subsonic stream bdry = %i!", bdry);



}


void RiemannInvariants::get_free_stream(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, double t_x, double t_y,
					double* new_variables, 
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, int& boundary, bool solid ){

	double lambda_1 = RiemannInvariants::get_ev1(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_2_3 = RiemannInvariants::get_ev2(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);
	double lambda_4 = RiemannInvariants::get_ev4(rho, rho_v_x, rho_v_y, rho_energy, n_x, n_y);


	double rho_new= 0.; double rho_v_x_new=0.; double rho_v_y_new=0.; double rho_energy_new=0.; 
if(solid==true){  //solid wall
			rho_new = rho;
			rho_energy_new = rho_energy;
			rho_v_x_new = rho_v_x - 2*n_x*( rho_v_x*n_x+ rho_v_y*n_y); 	
			rho_v_y_new = rho_v_y- 2*n_y*( rho_v_x*n_x+ rho_v_y*n_y); 
			boundary = 0;
	}else{
			if((lambda_1>=0)&&(lambda_2_3>=0)&&(lambda_4>=0)){ //supersonic outlet{

					rho_new= rho; 
					rho_v_x_new=rho_v_x; 	
					rho_v_y_new = rho_v_y;
					rho_energy_new = rho_energy;
					boundary = 1;
		}else if((lambda_1<0)&&(lambda_2_3<0)&&(lambda_4<0)){//supersonic inlet

					rho_new = rho_ext; 
					rho_v_x_new = rho_v_x_ext; 	
					rho_v_y_new = rho_v_y_ext;
					rho_energy_new = rho_energy_ext;
					boundary = 2;			
		}else{

	}

	new_variables[0]=  rho_new; 
	new_variables[1] = rho_v_x_new; 
	new_variables[2] = rho_v_y_new;
	new_variables[3] = rho_energy_new;
}
//-----------------------------------------------------
//-------------Discontinuity Detector---------------------
//-----------------------------------------------------

DiscontinuityDetector::DiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double>  > spaces, 
                                             Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : spaces(spaces), solutions(solutions)
{
  for(int i = 0; i < solutions.size(); i++)
    this->solutionsInternal.push_back((Solution<double>*)solutions[i].get());
};

DiscontinuityDetector::~DiscontinuityDetector()
{};

KrivodonovaDiscontinuityDetector::KrivodonovaDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double>  > spaces, 
                                                                   Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : DiscontinuityDetector(spaces, solutions)
{
  // A check that all meshes are the same in the spaces.
  unsigned int mesh0_seq = spaces[0]->get_mesh()->get_seq();
  for(unsigned int i = 0; i < spaces.size(); i++)
    if(spaces[i]->get_mesh()->get_seq() != mesh0_seq)
      throw Hermes::Exceptions::Exception("So far DiscontinuityDetector works only for single mesh.");
  mesh = spaces[0]->get_mesh();
};

KrivodonovaDiscontinuityDetector::~KrivodonovaDiscontinuityDetector()
{};

double KrivodonovaDiscontinuityDetector::calculate_h(Element* e, int polynomial_order)
{
  double h = std::sqrt(std::pow(e->vn[(0 + 1) % e->get_nvert()]->x - e->vn[0]->x, 2) + std::pow(e->vn[(0 + 1) % e->get_nvert()]->y - e->vn[0]->y, 2));
  for(int edge_i = 0; edge_i < e->get_nvert(); edge_i++) {
    double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->y - e->vn[edge_i]->y, 2));
    if(edge_length < h)
      h = edge_length;
  }
  return std::pow(h, (0.5 * (H2D_GET_H_ORDER(spaces[0]->get_element_order(e->id)) 
    + 
    H2D_GET_V_ORDER(spaces[0]->get_element_order(e->id)))
    + 1) / 2);
}

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids()
{
  return get_discontinuous_element_ids(1.0);
};

std::set<int>& KrivodonovaDiscontinuityDetector::get_discontinuous_element_ids(double threshold)
{
  Element* e;
if(!discontinuous_element_ids.empty()) discontinuous_element_ids.clear();
  for_all_active_elements(e, mesh)
  {
    bool element_inserted = false;
    for(int edge_i = 0; edge_i < e->get_nvert() && !element_inserted; edge_i++)
		{
      if(calculate_relative_flow_direction(e, edge_i) <= 0 && !e->en[edge_i]->bnd)
      { 
        double jumps[4];
        calculate_jumps(e, edge_i, jumps);
        double diameter_indicator = calculate_h(e, spaces[0]->get_element_order(e->id));
        double edge_length = std::sqrt(std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->x - e->vn[edge_i]->x, 2) + std::pow(e->vn[(edge_i + 1) % e->get_nvert()]->y - e->vn[edge_i]->y, 2));
        double norms[4];
        calculate_norms(e, edge_i, norms);

        // Number of component jumps tested.
        unsigned int component_checked_number = 1;
        for(unsigned int component_i = 0; component_i < component_checked_number; component_i++) {
          if(norms[component_i] < 1E-8)
            continue;
          double discontinuity_detector = jumps[component_i] / (diameter_indicator * edge_length * norms[component_i]);
          if(discontinuity_detector > threshold)
          {
            discontinuous_element_ids.insert(e->id);
            element_inserted = true;
            break;
          }
        }
      }
		}		
  }
  return discontinuous_element_ids;
};


double KrivodonovaDiscontinuityDetector::calculate_relative_flow_direction(Element* e, int edge_i)
{
  // Set active element to the two solutions (density_vel_x, density_vel_y).
  solutions[0]->set_active_element(e);
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[1]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 20, e->get_mode());
  double3* pt = solutions[1]->get_quad_2d()->get_points(eo, e->get_mode());
  int np = solutions[1]->get_quad_2d()->get_num_points(eo, e->get_mode());

  double3* tan;
  Geom<double>* geom = init_geom_surf(solutions[1]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);

  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2]* tan[i][2];

  // Calculate.
  Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
  Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);
  Func<double>* density				= init_fn(solutions[0].get(), eo);
	

  double result = 0.0;
  for(int point_i = 0; point_i < np; point_i++)
    result += jwt[point_i] * (density_vel_x->val[point_i]/density->val[point_i] * geom->nx[point_i] + density_vel_y->val[point_i]/density->val[point_i] * geom->ny[point_i]);

	
  geom->free();
  delete geom;
  delete [] jwt;
  density_vel_x->free_fn();
  density_vel_y->free_fn();
	density->free_fn();
	delete density;
  delete density_vel_x;
  delete density_vel_y;

  return result;
};

void KrivodonovaDiscontinuityDetector::calculate_jumps(Element* e, int edge_i, double result[4])
{
  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8, e->get_mode());
  double3* pt = solutions[0]->get_quad_2d()->get_points(eo, e->get_mode());
  int np = solutions[0]->get_quad_2d()->get_num_points(eo, e->get_mode());

  // Initialize the NeighborSearch.
  NeighborSearch<double> ns(e, mesh);
  ns.set_active_edge(edge_i);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Go through all neighbors.
  for(int neighbor_i = 0; neighbor_i < ns.get_num_neighbors(); neighbor_i++) {
    ns.set_active_segment(neighbor_i);

    // Set active element to the solutions.
    solutions[0]->set_active_element(e);
    solutions[1]->set_active_element(e);
    solutions[2]->set_active_element(e);
    solutions[3]->set_active_element(e);

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_central_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_central_transformations(neighbor_i, trf_i));
    }

    double3* tan;
    Geom<double>* geom = init_geom_surf(solutions[0]->get_refmap(), surf_pos.surf_num, surf_pos.marker, eo, tan);
    double* jwt = new double[np];
    for(int i = 0; i < np; i++)
      jwt[i] = pt[i][2] * tan[i][2];

    // Prepare functions on the central element.
    Func<double>* density = init_fn(solutions[0].get(), eo);
    Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
    Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);
    Func<double>* energy = init_fn(solutions[3].get(), eo);

    // Set neighbor element to the solutions.
    solutions[0]->set_active_element(ns.get_neighb_el());
    solutions[1]->set_active_element(ns.get_neighb_el());
    solutions[2]->set_active_element(ns.get_neighb_el());
    solutions[3]->set_active_element(ns.get_neighb_el());

    // Push all the necessary transformations.
    for(unsigned int trf_i = 0; trf_i < ns.get_neighbor_n_trans(neighbor_i); trf_i++) {
      solutions[0]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[1]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[2]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
      solutions[3]->push_transform(ns.get_neighbor_transformations(neighbor_i, trf_i));
    }

    // Prepare functions on the neighbor element.
    Func<double>* density_neighbor = init_fn(solutions[0].get(), eo);
    Func<double>* density_vel_x_neighbor = init_fn(solutions[1].get(), eo);
    Func<double>* density_vel_y_neighbor = init_fn(solutions[2].get(), eo);
    Func<double>* energy_neighbor = init_fn(solutions[3].get(), eo);

    DiscontinuousFunc<double> density_discontinuous(density, density_neighbor, true);
    DiscontinuousFunc<double> density_vel_x_discontinuous(density_vel_x, density_vel_x_neighbor, true);
    DiscontinuousFunc<double> density_vel_y_discontinuous(density_vel_y, density_vel_y_neighbor, true);
    DiscontinuousFunc<double> energy_discontinuous(energy, energy_neighbor, true);

    for(int point_i = 0; point_i < np; point_i++) {
      result[0] += jwt[point_i] * std::abs(density_discontinuous.val[point_i] - density_discontinuous.val_neighbor[point_i]); 
      result[1] += jwt[point_i] * std::abs(density_vel_x_discontinuous.val[point_i] - density_vel_x_discontinuous.val_neighbor[point_i]);
      result[2] += jwt[point_i] * std::abs(density_vel_y_discontinuous.val[point_i] - density_vel_y_discontinuous.val_neighbor[point_i]);
      result[3] += jwt[point_i] * std::abs(energy_discontinuous.val[point_i] - energy_discontinuous.val_neighbor[point_i]);

			/*result[0] += jwt[point_i] * (density_discontinuous.val[point_i] - density_discontinuous.val_neighbor[point_i]); 
      result[1] += jwt[point_i] * (density_vel_x_discontinuous.val[point_i] - density_vel_x_discontinuous.val_neighbor[point_i]);
      result[2] += jwt[point_i] * (density_vel_y_discontinuous.val[point_i] - density_vel_y_discontinuous.val_neighbor[point_i]);
      result[3] += jwt[point_i] * (energy_discontinuous.val[point_i] - energy_discontinuous.val_neighbor[point_i]);*/
    }

    geom->free();
    delete geom;
    delete [] jwt;
    density->free_fn();
    density_vel_x->free_fn();
    density_vel_y->free_fn();
    energy->free_fn();
    density_neighbor->free_fn();
    density_vel_x_neighbor->free_fn();
    density_vel_y_neighbor->free_fn();
    energy_neighbor->free_fn();

    delete density;
    delete density_vel_x;
    delete density_vel_y;
    delete energy;
    delete density_neighbor;
    delete density_vel_x_neighbor;
    delete density_vel_y_neighbor;
    delete energy_neighbor;
  }

  result[0] = std::abs(result[0]);
  result[1] = std::abs(result[1]);
  result[2] = std::abs(result[2]);
  result[3] = std::abs(result[3]);
};

void KrivodonovaDiscontinuityDetector::calculate_norms(Element* e, int edge_i, double result[4])
{
  // Set active element to the solutions.
  solutions[0]->set_active_element(e);
  solutions[1]->set_active_element(e);
  solutions[2]->set_active_element(e);
  solutions[3]->set_active_element(e);

  // The values to be returned.
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
  result[3] = 0.0;

  // Set Geometry.
  SurfPos surf_pos;
  surf_pos.marker = e->marker;
  surf_pos.surf_num = edge_i;

  int eo = solutions[0]->get_quad_2d()->get_edge_points(surf_pos.surf_num, 8, e->get_mode());
  int np = solutions[0]->get_quad_2d()->get_num_points(eo, e->get_mode());

  // Calculate.
  Func<double>* density = init_fn(solutions[0].get(), eo);
  Func<double>* density_vel_x = init_fn(solutions[1].get(), eo);
  Func<double>* density_vel_y = init_fn(solutions[2].get(), eo);
  Func<double>* energy = init_fn(solutions[3].get(), eo);

  for(int point_i = 0; point_i < np; point_i++) {
    result[0] = std::max(result[0], std::abs(density->val[point_i]));
    result[1] = std::max(result[1], std::abs(density_vel_x->val[point_i]));
    result[2] = std::max(result[2], std::abs(density_vel_y->val[point_i]));
    result[3] = std::max(result[3], std::abs(energy->val[point_i]));
  }



  density->free_fn();
  density_vel_x->free_fn();
  density_vel_y->free_fn();
  energy->free_fn();

  delete density;
  delete density_vel_x;
  delete density_vel_y;
  delete energy;
};



