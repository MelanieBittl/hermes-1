class CustomInitialCondition_rho : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_rho(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
  ~CustomInitialCondition_rho(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition_rho(this->mesh);
  }
};

class CustomInitialCondition_v_x : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_v_x(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
  ~CustomInitialCondition_v_x(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition_v_x(this->mesh);
  }
};

class CustomInitialCondition_v_y : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_v_y(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};
  ~CustomInitialCondition_v_y(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  virtual Ord ord(double x, double y) const ;

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition_v_y(this->mesh);
  }
};

class CustomInitialCondition_e : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition_e(MeshSharedPtr mesh, double kappa) : ExactSolutionScalar<double>(mesh), kappa(kappa) {};
  ~CustomInitialCondition_e(){};

  virtual void derivatives (double x, double y, double& dx, double& dy) const ;

  virtual double value (double x, double y) const;

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition_e(this->mesh, this->kappa);
  }

  virtual Ord ord(double x, double y) const ;
  double kappa;
};

void CustomInitialCondition_rho::derivatives(double x, double y, double& dx, double& dy) const {      
  dx = 0.0;
  dy = 0.0;
};

double CustomInitialCondition_rho::value(double x, double y) const {       
  if(x< (1./6.+y/std::sqrt(3.))) return 8.0;
  else			return 1.4;

};

Ord CustomInitialCondition_rho::ord(double x, double y) const {
  return Ord(2);
};

void CustomInitialCondition_v_x::derivatives(double x, double y, double& dx, double& dy) const {      
  dx = 0.0;
  dy = 0.0;
};

double CustomInitialCondition_v_x::value(double x, double y) const {       
  if(x< (1./6.+y/std::sqrt(3.))) return 8.25*std::cos(M_PI/6.)*8.0;
  else			return 0.0;

};

Ord CustomInitialCondition_v_x::ord(double x, double y) const {
  return Ord(2);
};
void CustomInitialCondition_v_y::derivatives(double x, double y, double& dx, double& dy) const {      
  dx = 0.0;
  dy = 0.0;
};

double CustomInitialCondition_v_y::value(double x, double y) const {       
  if(x< (1./6.+y/std::sqrt(3.))) return -8.25*std::sin(M_PI/6.)*8.0;
  else			return 0.0;

};

Ord CustomInitialCondition_v_y::ord(double x, double y) const {
  return Ord(2);
};



void CustomInitialCondition_e::derivatives(double x, double y, double& dx, double& dy) const {      
  dx = 0.0;
  dy = 0.0;
};

double CustomInitialCondition_e::value(double x, double y) const {       
  if(x< (1./6.+y/std::sqrt(3.))) return QuantityCalculator::calc_energy(8.0, 8.25*std::cos(M_PI/6.)*8.0 ,-8.25*std::sin(M_PI/6.)*8.0, 116.5, kappa);
  else			return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);
};

Ord CustomInitialCondition_e::ord(double x, double y) const {
  return Ord(2);
};