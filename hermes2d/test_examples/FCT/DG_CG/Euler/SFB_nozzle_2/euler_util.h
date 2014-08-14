#ifndef EULER_UTIL_H
#define EULER_UTIL_H

#include "hermes2d.h"





using namespace Hermes;
using namespace Hermes::Hermes2D;





//------------------------- Class calculating various quantities-----------
class QuantityCalculator
{
public:
  // Calculates energy*rho from other quantities. 
  static double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure, double kappa);
 
  // Calculates pressure from other quantities.
  static double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
 
  // Calculates speed of sound.
  static double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);

  // Calculates enthalpy.
  static double enthalpy(double rho, double rho_v_x, double rho_v_y, double energy, double kappa);
};

//----------------------------------------------------------------------------------
//----------------------Filters.--------------------------------------
//----------------------------------------------------------------------------------
class MachNumberFilter : public Hermes::Hermes2D::SimpleFilter<double>
{
public: 
  MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
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
    MachNumberFilter* filter = new MachNumberFilter(slns, this->kappa);
          
    return filter;
  }

protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
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
  PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa) : SimpleFilter<double>(solutions), kappa(kappa) {};
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
    PressureFilter* filter = new PressureFilter(slns, this->kappa);

    return filter;

  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa;
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
  EntropyFilter(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, double kappa, double rho_ext, double p_ext) : SimpleFilter<double>(solutions), kappa(kappa), rho_ext(rho_ext), p_ext(p_ext) {};
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
    EntropyFilter* filter = new EntropyFilter(slns, this->kappa, rho_ext, p_ext);

    return filter;
  }
protected:
  virtual void filter_fn(int n, Hermes::vector<double*> values, double* result);

  double kappa, rho_ext, p_ext;
};

//----------------------------------------------------------------------------------
//----------------RiemannInvariants---------------------
//----------------------------------------------------------------------------------
class RiemannInvariants
{
public:
	RiemannInvariants(double kappa) : kappa(kappa){};
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

	void get_du_du(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, int bdry, int entry_j, double* dudu);

	void get_ghost_state(int bdry,double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y,double t_x, double t_y,
					double rho_ext, double rho_v_x_ext, double rho_v_y_ext, double rho_energy_ext, double* ghost_state, bool solid);

protected:
	double kappa;

};

//----------------------------------------------------------------------------------
//-----------------------DiscontinuityDetector------------------------------
//----------------------------------------------------------------------------------
class DiscontinuityDetector
{
public:
  /// Constructor.
  DiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, 
    Hermes::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
  ~DiscontinuityDetector();

  /// Return a reference to the inner structures.
  virtual std::set<int>& get_discontinuous_element_ids() = 0;
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRho() { return this->oscillatory_element_idsRho; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoVX() { return this->oscillatory_element_idsRhoVX; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoVY() { return this->oscillatory_element_idsRhoVY; }
  std::set<std::pair<int, double> >& get_oscillatory_element_idsRhoE() { return this->oscillatory_element_idsRhoE; }

protected:
  /// Members.
  Hermes::vector<SpaceSharedPtr<double> > spaces;
  Hermes::vector<MeshFunctionSharedPtr<double> > solutions;
  Hermes::vector<Solution<double>*> solutionsInternal;
  std::set<int> discontinuous_element_ids;
  std::set<std::pair<int, double> > oscillatory_element_idsRho;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVX;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoVY;
  std::set<std::pair<int, double> > oscillatory_element_idsRhoE;
  MeshSharedPtr mesh;
};

class KrivodonovaDiscontinuityDetector : public DiscontinuityDetector
{
public:
  /// Constructor.
  KrivodonovaDiscontinuityDetector(Hermes::vector<SpaceSharedPtr<double> > spaces, 
    Hermes::vector<MeshFunctionSharedPtr<double> > solutions);

  /// Destructor.
  ~KrivodonovaDiscontinuityDetector();

  /// Return a reference to the inner structures.
  std::set<int>& get_discontinuous_element_ids();
  std::set<int>& get_discontinuous_element_ids(double threshold);



protected:
  /// Calculates relative (w.r.t. the boundary edge_i of the Element e).
  double calculate_relative_flow_direction(Element* e, int edge_i);

  /// Calculates jumps of all solution components across the edge edge_i of the Element e.
  void calculate_jumps(Element* e, int edge_i, double result[4]);

  /// Calculates h.
  double calculate_h(Element* e, int polynomial_order);

  /// Calculates the norm of the solution on the central element.
  void calculate_norms(Element* e, int edge_i, double result[4]);

};

//----------------------------------------------------------------------------------------
///-------------------------------------------- for boundary calculations
////-----------

class Boundary_helpers
{
	public:

		static double calculate_A_n_U(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
							double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry){

			double u[4] = {rho_new-rho, rho_v_x_new-rho_v_x, rho_v_y_new-rho_v_y, rho_energy_new-rho_energy };

			if((u[0]==0.)&&(u[1]==0.)&&(u[2]==0.)&&(u[3]==0.)) return 0.;

			double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
			double v_n = v_x_mean*n_x+v_y_mean*n_y;
			double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
			double b_2 = (kappa-1)/(c_mean*c_mean);
			double b_1 = b_2*q;
			double lambda[4] ={fabs(v_n- c_mean),fabs(v_n),fabs(v_n + c_mean),fabs(v_n)};

		//Eintraege siehe Diss Moeller Appendix C
			double R[4][4] = { 1, 1, 1, 0, 
												v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
												v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
												H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

			double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
												1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
												0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
												n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

		/*u[0] = -rho;
		u[1]= -rho_v_x;
		u[2]= -rho_v_y;
			u[3] = -rho_energy ;*/

			double result =0.;

			for(int i =0;i<4;i++)
					for(int j=0;j<4;j++)
							result +=R[entry][j]*L[j][i]*u[i]*lambda[j];


			return result;		


		};





		static void calculate_A_n(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
							double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry_i, double* A_n){
			double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));


			double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
			double v_n = v_x_mean*n_x+v_y_mean*n_y;
			double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
			double b_2 = (kappa-1)/(c_mean*c_mean);
			double b_1 = b_2*q;
			double lambda[4] ={fabs(v_n- c_mean),fabs(v_n),fabs(v_n + c_mean),fabs(v_n)};



		//Eintraege siehe Diss Moeller Appendix C
			double R[4][4] = { 1, 1, 1, 0, 
												v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
												v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
												H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

			double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
												1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
												0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
												n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

			for(int i =0;i<4;i++) A_n[i]=0.;
	
			for(int j=0;j<4;j++){
					A_n[0] +=R[entry_i][j]*L[j][0]*lambda[j];
					A_n[1] +=R[entry_i][j]*L[j][1]*lambda[j];
					A_n[2] +=R[entry_i][j]*L[j][2]*lambda[j];
					A_n[3] +=R[entry_i][j]*L[j][3]*lambda[j];
				}



		};






		static void calculate_P_plus(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
							double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry_i, double* A_n){
			//double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			//double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			//double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));

double v_x_mean = (rho_v_x/rho + rho_v_x_new/rho_new)*0.5;
double v_y_mean = (rho_v_y/rho + rho_v_y_new/rho_new)*0.5;
double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)+QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa))*0.5;

			double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
			double v_n = v_x_mean*n_x+v_y_mean*n_y;
			double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
			double b_2 = (kappa-1)/(c_mean*c_mean);
			double b_1 = b_2*q;

			double lambda[4] ={ ((v_n- c_mean)>0? v_n- c_mean:0. ) ,(v_n>0? v_n: 0. ),((v_n+ c_mean)>0? v_n+ c_mean:0. ),(v_n>0? v_n: 0. )};



		//Eintraege siehe Diss Moeller Appendix C
			double R[4][4] = { 1, 1, 1, 0, 
												v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
												v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
												H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

			double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
												1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
												0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
												n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

			for(int i =0;i<4;i++) A_n[i]=0.;
	
			for(int j=0;j<4;j++){
					A_n[0] +=R[entry_i][j]*L[j][0]*lambda[j];
					A_n[1] +=R[entry_i][j]*L[j][1]*lambda[j];
					A_n[2] +=R[entry_i][j]*L[j][2]*lambda[j];
					A_n[3] +=R[entry_i][j]*L[j][3]*lambda[j];
				}



		};

	static void calculate_P_minus(double rho, double rho_v_x, double rho_v_y, double rho_energy, double n_x, double n_y, 
							double rho_new, double rho_v_x_new, double rho_v_y_new, double rho_energy_new, double kappa, int entry_i, double* A_n){
			//double v_x_mean = (rho_v_x/std::sqrt(rho) + rho_v_x_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			//double v_y_mean = (rho_v_y/std::sqrt(rho) + rho_v_y_new/std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));
			//double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)*std::sqrt(rho) + QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa)*std::sqrt(rho_new))/(std::sqrt(rho) +std::sqrt(rho_new));


double v_x_mean = (rho_v_x/rho + rho_v_x_new/rho_new)*0.5;
double v_y_mean = (rho_v_y/rho + rho_v_y_new/rho_new)*0.5;
double H_mean = (QuantityCalculator::enthalpy(rho, rho_v_x,rho_v_y, rho_energy, kappa)+QuantityCalculator::enthalpy(rho_new, rho_v_x_new,rho_v_y_new, rho_energy_new, kappa))*0.5;


			double c_mean = std::sqrt((kappa-1)*(H_mean- 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean)));
			double v_n = v_x_mean*n_x+v_y_mean*n_y;
			double q = 0.5*(v_x_mean*v_x_mean+v_y_mean*v_y_mean);
			double b_2 = (kappa-1)/(c_mean*c_mean);
			double b_1 = b_2*q;

			double lambda[4] ={ ((v_n- c_mean)<0? v_n- c_mean:0. ) ,(v_n<0? v_n: 0. ),((v_n+ c_mean)<0? v_n+ c_mean:0. ),(v_n<0? v_n: 0. )};



		//Eintraege siehe Diss Moeller Appendix C
			double R[4][4] = { 1, 1, 1, 0, 
												v_x_mean-c_mean*n_x, v_x_mean, v_x_mean+c_mean*n_x, n_y,
												v_y_mean-c_mean*n_y, v_y_mean, v_y_mean+c_mean*n_y, -n_x,
												H_mean-c_mean*v_n, q, H_mean+c_mean*v_n, v_x_mean*n_y-v_y_mean*n_x};

			double L[4][4] = {0.5*(b_1+v_n/c_mean), 0.5*(-b_2*v_x_mean-n_x/c_mean), 0.5*(-b_2*v_y_mean-n_y/c_mean), 0.5*b_2,
												1-b_1, b_2*v_x_mean, b_2*v_y_mean, -b_2,
												0.5*(b_1-v_n/c_mean), 0.5*(-b_2*v_x_mean+n_x/c_mean), 0.5*(-b_2*v_y_mean+n_y/c_mean), 0.5*b_2,
												n_x*v_y_mean-n_y*v_x_mean, n_y,  -n_x, 0};

			for(int i =0;i<4;i++) A_n[i]=0.;
	
			for(int j=0;j<4;j++){
					A_n[0] +=R[entry_i][j]*L[j][0]*lambda[j];
					A_n[1] +=R[entry_i][j]*L[j][1]*lambda[j];
					A_n[2] +=R[entry_i][j]*L[j][2]*lambda[j];
					A_n[3] +=R[entry_i][j]*L[j][3]*lambda[j];
				}

		};



};


#endif
