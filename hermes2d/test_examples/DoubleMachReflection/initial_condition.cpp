class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh, int component, double kappa) : ExactSolutionScalar<double>(mesh), component(component), kappa(kappa)
  {
    sqrt3 = std::sqrt(3.0);
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const 
  {
    dx = 0.0;
    dy = 0.0;
  }

  virtual double value (double x, double y) const
  {
    switch(this->component)
    {
    case 0:
      if(x< (1./6. + ((y * (1 + 20 * this->time) / sqrt3))))
        return 8.0;
      else		
        return 1.4;
      break;
    case 1:
      if(x< (1./6. + ((y * (1 + 20 * this->time) / sqrt3))))
        return 8.25*std::cos(M_PI/6.)*8.0;
      else		
        return 0.0;
      break;
    case 2:
      if(x< (1./6. + ((y * (1 + 20 * this->time) / sqrt3))))
        return -8.25*std::sin(M_PI/6.)*8.0;
      else		
        return 0.0;
      break;
    case 3:
      if(x< (1./6. + ((y * (1 + 20 * this->time) / sqrt3))))
        return QuantityCalculator::calc_energy(8.0, 8.25*std::cos(M_PI/6.)*8.0 ,-8.25*std::sin(M_PI/6.)*8.0, 116.5, kappa);
      else		
        return QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);
      break;
    }
  }

  virtual Ord ord(double x, double y) const
  {
    return Ord(2);
  }

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition(this->mesh, this->component, this->kappa);
  }

  int component;
  double time;
  double sqrt3;
  double kappa;
};