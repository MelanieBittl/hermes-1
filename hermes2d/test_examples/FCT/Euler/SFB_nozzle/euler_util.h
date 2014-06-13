#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

//------------------------- Class calculating various quantities-----------
class QuantityCalculator
{
public:
  // Calculates energy from other quantities.
  static double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double GAMMA);
 
  // Calculates pressure from other quantities.
  static double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double GAMMA);
 
  // Calculates speed of sound.
  static double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double GAMMA);

  // Calculates enthalpy.
  static double enthalpy(double rho, double rho_v_x, double rho_v_y, double energy, double GAMMA);
 // Calculates mach number.
  static double calc_mach(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
};


//----------------------Filters.--------------------------------------

class MachNumberFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double gamma) : SimpleFilter<double>(solutions), gamma(gamma) {};
  ~MachNumberFilter() 
  {
  };

   MeshFunction<double>*  clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    MachNumberFilter* filter = new MachNumberFilter(slns, this->gamma);
          
    return filter;
  }

protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double gamma;
};



class VelocityFilter_x : public Hermes::Hermes2D::SimpleFilter<double>
{
public:
  VelocityFilter_x(Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : SimpleFilter<double>(solutions) {};
  ~VelocityFilter_x() 
  {
  };

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    VelocityFilter_x* filter = new VelocityFilter_x(slns);
    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
};


class VelocityFilter_y : public Hermes::Hermes2D::SimpleFilter<double>
{
public:

  VelocityFilter_y(Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : SimpleFilter<double>(solutions) {};
  ~VelocityFilter_y() 
  {
  };

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }

    VelocityFilter_y* filter = new VelocityFilter_y(slns);
    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);
};

class PressureFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double gamma) : SimpleFilter<double>(solutions), gamma(gamma) {};
  ~PressureFilter() 
  {
  };

  MeshFunction<double>*  clone() const
  {	
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    PressureFilter* filter = new PressureFilter(slns, this->gamma);

    return filter;

  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double gamma;
};
class RadiusVelocityFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  RadiusVelocityFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions) : SimpleFilter<double>(solutions) {};

  MeshFunction<double>* clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    RadiusVelocityFilter* filter = new RadiusVelocityFilter(slns);
    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

};

class EntropyFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  EntropyFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double gamma, double rho_ext, double p_ext) : SimpleFilter<double>(solutions), gamma(gamma), rho_ext(rho_ext), p_ext(p_ext) {};
  ~EntropyFilter() 
  {
  };
  MeshFunction<double>*  clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
    EntropyFilter* filter = new EntropyFilter(slns, this->gamma, rho_ext, p_ext);

    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double gamma, rho_ext, p_ext;
};


class TempFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  TempFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double R, double gamma) : SimpleFilter<double>(solutions), R(R),gamma(gamma) {};
  ~TempFilter() 
  {
  };

   MeshFunction<double>*  clone() const
  {
    Hermes::vector<MeshFunctionSharedPtr<double> > slns;
    for(int i = 0; i < this->num; i++)
    {
      Solution<double>* solution = dynamic_cast<Solution<double>*>(this->sln[i].get());
      if(solution && solution->get_type() == HERMES_SLN)
      {
        slns.push_back(new Solution<double>());
        slns.back()->copy(this->sln[i]);
      }
      else
        slns.push_back(this->sln[i]->clone()); 
    }
   TempFilter* filter = new TempFilter(slns, this->R, this->gamma);
          
    return filter;
  } 
  protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double R, gamma;
};

//----------------RiemannInvariants---------------------
class RiemannInvariants
{
public:
	RiemannInvariants(double GAMMA) : GAMMA(GAMMA){};
	double get_w1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);
	double get_w2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);
	double get_w3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double t_x, double t_y);
	double get_w4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);

	bool solid_wall(double rho, double rho_v_x, double rho_v_y, double n_x, double n_y);

	double get_ev1(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);
	double get_ev2(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);
	double get_ev3(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);
	double get_ev4(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y);

	double get_pressure(double w_1, double w_2, double w_3, double w_4);
	double get_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x);
	double get_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y);
	double get_speed_sound(double w_1, double w_4);

	double get_rho(double w_1, double w_2, double w_3, double w_4);
	double get_rho_v_x(double w_1, double w_2, double w_3, double w_4, double n_x, double t_x);
	double get_rho_v_y(double w_1, double w_2, double w_3, double w_4,double n_y, double t_y);
	double get_energy(double w_1, double w_2, double w_3, double w_4, double n_x, double n_y, double t_x, double t_y);

	void get_free_stream(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, double t_x, double t_y,
					double* new_variables, 
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, int& boundary, bool solid );

	int get_bdry_info(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid);

void get_ghost_state(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state);

	void get_du_du(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, int bdry, int entry_j, double* dudu);

protected:
	double GAMMA;

};


//--------------------Euler-Fluxes-------------------------
class EulerFluxes
{
public:
  EulerFluxes(double GAMMA) : GAMMA(GAMMA) {}

  template<typename Scalar>
  Scalar A_1_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {

    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_1_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_2_0_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(1.0);
  }

  template<typename Scalar>
  Scalar A_2_0_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(0.0);
  }

  template<typename Scalar>
  Scalar A_1_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy) {
    return Scalar(- ((rho_v_x * rho_v_x) / (rho * rho)) + 0.5 * (GAMMA - 1.0) * 
            ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_1_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3. - GAMMA) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - GAMMA) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_1_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(GAMMA - 1.);
  }

  template<typename Scalar>
  Scalar A_2_1_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_1_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_2_1_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_1_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- rho_v_x * rho_v_y / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_y / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(rho_v_x / rho);
  }

  template<typename Scalar>
  Scalar A_1_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(0);
  }

  template<typename Scalar>
  Scalar A_2_2_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(- ((rho_v_y * rho_v_y) / (rho * rho)) + 0.5 * (GAMMA - 1.0) 
            * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) /   (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_2_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - GAMMA) * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((3.0 - GAMMA) * (rho_v_y / rho));
  }

  template<typename Scalar>
  Scalar A_2_2_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(GAMMA - 1.);
  }

  template<typename Scalar>
  Scalar A_1_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((rho_v_x / rho) * (((GAMMA - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho)))
      - (GAMMA * energy / rho)));
  }

  template<typename Scalar>
  Scalar A_1_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((GAMMA * energy / rho) - (GAMMA - 1.0) * rho_v_x * rho_v_x / (rho * rho)
      - 0.5 * (GAMMA - 1.0) * (rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - GAMMA) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_1_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(GAMMA * (rho_v_x / rho));
  }

  template<typename Scalar>
  Scalar A_2_3_0(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    /*return Scalar(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (GAMMA - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho))) + (rho_v_y / rho) 
            * (GAMMA - 1.0) * ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (2 * rho * rho)));*/

		return Scalar(- (rho_v_y * energy) / (rho * rho) - (rho_v_y / (rho * rho)) * (GAMMA - 1.0) 
            * (energy - ((rho_v_x * rho_v_x + rho_v_y * rho_v_y) / (rho))));


  }

  template<typename Scalar>
  Scalar A_2_3_1(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((1.0 - GAMMA) * (rho_v_x * rho_v_y) / (rho * rho));
  }

  template<typename Scalar>
  Scalar A_2_3_2(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar((energy / rho) + (1 / rho) * (GAMMA - 1.0) * ( energy - ((rho_v_x * rho_v_x 
            + rho_v_y * rho_v_y) / (2 * rho))) + (1.0 - GAMMA) * ((rho_v_y * rho_v_y) / (rho * rho)));
  }

  template<typename Scalar>
  Scalar A_2_3_3(Scalar rho, Scalar rho_v_x, Scalar rho_v_y, Scalar energy){
    return Scalar(GAMMA * rho_v_y / rho);
  }
  protected:
    double GAMMA;
};

#endif
