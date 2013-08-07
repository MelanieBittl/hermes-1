#ifndef __SlopeLimiter_Solution_H
#define __SlopeLimiter_Solution_H
#include "hermes2d.h" 
 #include "shapeset_taylor.h"


using namespace Hermes;
using namespace Hermes::Hermes2D;


 class KuzminOscillationDetector
{
public:
  /// Constructor.
  KuzminOscillationDetector(SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double> sln_fct);

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
    MeshSharedPtr mesh;
  SpaceSharedPtr<double>  space;
  MeshFunctionSharedPtr<double> solution;
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
 	SlopeLimiterSolution(MeshFunctionSharedPtr<double> sln, SpaceSharedPtr<double> space);
 	~SlopeLimiterSolution(){ delete detector;};
 	
 	void limit_solution_according_to_detector(bool p1_only=false);
 	
 	double* get_alpha_first_order(){return detector->get_alpha_first_order();};
 	double* get_alpha_second_order(){return detector->get_alpha_second_order();};
 	
 	void set_diff_linear();

 	
 protected:   
   MeshFunctionSharedPtr<double>  ref_sln;
   SpaceSharedPtr<double> space;
   KuzminOscillationDetector* detector;
 
 };












#endif
