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
   

  void compute_all(bool only_p1 = false);
  double* get_alpha_first_order(){return alpha_i_first_order;};
  double* get_alpha_second_order(){return alpha_i_second_order;};  

  void get_delta(Element* e, double &delta_x, double &delta_y); 

  
    /// Center.
  void find_centroid_values(Element* e, double &u_c);
  void find_centroid_derivatives(Element* e, double &u_dx_c, double &u_dy_c);
  void find_second_centroid_derivatives(Element* e, double &u_dxx_c, double &u_dxy_c, double &u_dyy_c);

protected:
  
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
 	~SlopeLimiterSolution(){ delete detector;};
 	
 	void limit_solution_according_to_detector(bool p1_only=false);
 	
 	double* get_alpha_first_order(){return detector->get_alpha_first_order();};
 	double* get_alpha_second_order(){return detector->get_alpha_second_order();};
 	
 	void set_diff_linear();
 	
 protected:
   
   Solution<double>* ref_sln;
   Space<double>* space;
   KuzminOscillationDetector* detector;
 
 };

 class RegSolution : public Solution<double>
 { 
 public:
 	RegSolution(Space<double>* space);
 	void linear_approx(double* u_c, Solution<double>* R_h_1, Solution<double>* R_h_2);
 	void linear_approx_d(double* u_c_d, Solution<double>* R_h); 	
 protected:   
   Space<double>* space;
 
 };











#endif
