#ifndef __SlopeLimiter_Solution_H
#define __SlopeLimiter_Solution_H
#include "hermes2d.h" 
#include "util.h"
 #include "shapeset_taylor.h"
 #include "l2_semi_cg_space.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;


 class KuzminOscillationDetector
{
public:
  /// Constructor.
  KuzminOscillationDetector(Space<double> * space, Solution<double> * solution);

  /// Destructor.
   ~KuzminOscillationDetector();
   

  void compute_all();
  double* get_alpha_first_order(){return alpha_i_first_order;};
  double* get_alpha_second_order(){return alpha_i_second_order;};
  
  double get_u_c(int id) {return u_c[id];}; 
  double get_u_dx_K(int id) {return u_dx_K[id];}; 
  double get_u_dy_K(int id){return u_dy_K[id];};  
  double get_u_dxx_K(int id){return u_dxx_K[id];}; 
  double get_u_dyy_K(int id){return u_dyy_K[id];};   
  double get_u_dxy_K(int id){return u_dxy_K[id];};   
  
  double get_delta_x(Element* e);
  double get_delta_y(Element* e);
  
  double get_alpha_first_order(Element* e);
  double get_alpha_second_order(Element* e);
  void get_alpha_all_orders(Element* e, double* alpha,double &u_c, double &u_dx_c, double &u_dy_c, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c);
  
    /// Center.
  void find_centroid_values(Element* e, double &u_c);
  void find_centroid_derivatives(Element* e, double &u_dx_c, double &u_dy_c);
  void find_second_centroid_derivatives(Element* e, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c);

protected:

  /// Vertices.
  void find_vertex_values(Element* e, double vertex_values[4]);
  void find_vertex_derivatives(Element* e, double vertex_derivatives[4][2]);  

	double find_value(Element* e, double x, double y, int derivative);

  /// Logic - 1st order.
  void find_u_i_min_max_first_order(Element* e, double u_i_min[4], double u_i_max[4]);
  void find_alpha_i_first_order(double u_i_min[4], double u_i_max[4], double u_c, double u_i[4], double &alpha_i);
  


  /// Logic - 2nd order.
  void find_u_i_min_max_second_order(Element* e, double u_d_i_min[4][2], double u_d_i_max[4][2]);
  void find_alpha_i_second_order(double u_d_i_min[4][2], double u_d_i_max[4][2], double &u_dx_c, double &u_dy_c, double u_d_i[4][2], double &alpha_i);
  
  
  /// Members.
    Mesh* mesh;
  Space<double> * space;
  Solution<double> * solution;
    VolumeAverage* vol_average;       

  
  double* u_c; 
  double* u_dx_K; 
  double* u_dy_K;  
  double* u_dxx_K; 
  double* u_dyy_K;   
  double* u_dxy_K;   
  
  double* alpha_i_first_order; 
  double* alpha_i_second_order;
  

  
  double* u_min_i;
  double* u_max_i;
  double* u_dx_min_i;
  double* u_dx_max_i;
  double* u_dy_min_i;
  double* u_dy_max_i;

    
};





  class Quad2DCheb;

 class SlopeLimiterSolution : public Solution<double>
 { 
 public:
 	SlopeLimiterSolution(Solution<double>* ref_sln, Space<double>* space);
 	
 	void limit_solution_according_to_detector();
 	
 protected:
   
   Solution<double>* ref_sln;
   Space<double>* space;
   KuzminOscillationDetector* detector;
 
 };













#endif
