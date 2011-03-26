#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

using namespace Laplace;

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  // Problem parameters.
  CustomWeakFormPoissonNewton(double h, double T0, std::string natural_bc_bnd_part) : WeakForm(1)
  {
    add_matrix_form(new VolumetricMatrixForms::DefaultMatrixFormStiffness(0, 0));
    add_matrix_form_surf(new SurfaceMatrixForms::DefaultSurfaceMatrixForm(0, 0, natural_bc_bnd_part, h));
    add_vector_form_surf(new SurfaceVectorForms::DefaultSurfaceVectorForm(0, natural_bc_bnd_part, h * T0));
  };
};

