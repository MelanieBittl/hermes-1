#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const bool BASE_VISUALIZATION = true;              // Set to "true" to enable base functions output.
const int P_INIT = 6;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 6;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double VOLUME_HEAT_SRC = 1.0;          // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 10.0;					// Fixed temperature on the boundary.

const Hermes::Hermes2D::ElementMode2D elementMode = HERMES_MODE_TRIANGLE;

double matrixFunction(int np, double* wt, Func<double>* u, Func<double>* v)
{
	double result = 0;
  for (int i = 0; i < np; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	return result;
}

double rhsValue = VOLUME_HEAT_SRC;

double rhsFunction(int np, double* wt, Func<double>* v)
{
	double result = 0;
  for (int i = 0; i < np; i++)
    result += wt[i] * v->val[i] * rhsValue;
	return result;
}

// Globals.
Hermes::Hermes2D::H1Space<double>* space;
int max_index;
Hermes::Hermes2D::Mesh* mesh;

// Storage.
double** A_values;
double* rhs_values;

// Cache loading.
void loadCache();

// Problem loading.
void loadProblemData();

// Cache calculation.
void calculateCache(CustomWeakFormPoisson& wf);

// Calculation from cache.
void calculateResultFromCache(CustomWeakFormPoisson& wf);
double handleDirichlet(CustomWeakFormPoisson& wf, int shape_indexBasis, int shape_indexDirichlet, Element* e);

// For comparison: calculation by using standard assembling.
void calculateResultAssembling(CustomWeakFormPoisson& wf);

int main(int argc, char* argv[])
{
  // Initialize the weak formulation.
  CustomWeakFormPoisson wf(new Hermes::Hermes2DFunction<double>(VOLUME_HEAT_SRC));

  Hermes2DApi.set_integral_param_value(numThreads, 1);

	//calculateCache(wf);

  loadProblemData();

	A_values = new double*[max_index];
  for(unsigned int i = 0; i < max_index; i++)
    A_values[i] = new double[max_index];

  rhs_values = new double[max_index];

  //loadCache();

  //calculateResultFromCache(wf);
	calculateResultAssembling(wf);
}

void loadCache()
{
	std::stringstream ssMatrix;
	std::stringstream ssRhs;
	ssMatrix << "matrixCache";
	ssRhs << "rhsCache";
	if(elementMode == HERMES_MODE_TRIANGLE)
	{
		ssMatrix << "Triangle";
		ssRhs << "Triangle";
	}
	else
	{
		ssMatrix << "Quad";
		ssRhs << "Quad";
	}

	std::ifstream matrixFormIn(ssMatrix.str());
	std::ifstream rhsFormIn(ssRhs.str());

	int index_i, index_j;
	int counter = 0;
  double valueTemp;
  while(matrixFormIn.good())
  {
    matrixFormIn >> index_i >> index_j >> valueTemp;
    A_values[index_i][index_j] = valueTemp;
    counter++;
  }

  counter = 0;
  while(rhsFormIn.good())
  {
    rhsFormIn >> index_i >> valueTemp;
    rhs_values[index_i] = valueTemp;
    counter++;
  }
}

void loadProblemData()
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mesh = new Mesh();
	std::stringstream ss;
	ss << "domain";
	if(elementMode == HERMES_MODE_TRIANGLE)
		ss << "Triangle";
	else
		ss << "Quad";
	ss << ".xml";
	try
	{
		mloader.load(ss.str().c_str(), mesh);
	}
	catch(std::exception& e)
	{
		std::cout << e.what();
	}

	for(int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Create a H1 space with default shapeset.
  Hermes::Hermes2D::DefaultEssentialBCConst<double>* bc_essential = new Hermes::Hermes2D::DefaultEssentialBCConst<double>("Bnd", FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double>* bcs = new Hermes::Hermes2D::EssentialBCs<double>(bc_essential);

  // Create an H1 space with default shapeset.
  space = new Hermes::Hermes2D::H1Space<double>(mesh, bcs, P_INIT);
}

void calculateCache(CustomWeakFormPoisson& wf)
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mesh = new Mesh();
	std::stringstream ss;
	ss << "domainCacheComputation";
	if(elementMode == HERMES_MODE_TRIANGLE)
		ss << "Triangle";
	else
		ss << "Quad";
	ss << ".xml";

	try
	{
		mloader.load(ss.str().c_str(), mesh);
	}
	catch(std::exception& e)
	{
		std::cout << e.what();
	}

	std::stringstream ssMatrix;
	std::stringstream ssRhs;
	ssMatrix << "matrixCache";
	ssRhs << "rhsCache";

	if(elementMode == HERMES_MODE_TRIANGLE)
	{
		ssMatrix << "Triangle";
		ssRhs << "Triangle";
	}
	else
	{
		ssMatrix << "Quad";
		ssRhs << "Quad";
	}

	std::ofstream matrixFormOut(ssMatrix.str());
	std::ofstream rhsFormOut(ssRhs.str());

  H1Shapeset shapeset;
  max_index = shapeset.get_max_index(elementMode) + 1;
	PrecalcShapeset uP(&shapeset);
	PrecalcShapeset vP(&shapeset);
	uP.set_active_element(mesh->get_element(0));
	vP.set_active_element(mesh->get_element(0));
  
	int np = g_quad_2d_std.get_num_points(elementMode == HERMES_MODE_TRIANGLE ? 20 : 24, elementMode);
	double3* pt = g_quad_2d_std.get_points(elementMode == HERMES_MODE_TRIANGLE ? 20 : 24, elementMode);
	double* jwt = new double[np];
	for(int i = 0; i < np; i++)
		jwt[i] = pt[i][2];

	RefMap rm;
	rm.set_active_element(mesh->get_element(0));
	
	for(int i = 0; i < max_index; i++)
  {
    if(i > 0)
      rhsFormOut << std::endl;
		vP.set_active_shape(i);
		Func<double>* v = init_fn(&vP, &rm, (elementMode == HERMES_MODE_TRIANGLE ? 20 : 24));

    for(int j = 0; j < max_index; j++)
    {
      if(i > 0 || j > 0)
        matrixFormOut << std::endl;
			uP.set_active_shape(j);
			Func<double>* u = init_fn(&uP, &rm, (elementMode == HERMES_MODE_TRIANGLE ? 20 : 24));
      matrixFormOut << i << ' ' << j << ' ' << matrixFunction(np, jwt, u, v);
    }
    rhsFormOut << i << ' ' << rhsFunction(np, jwt, v);
  }

  matrixFormOut.close();
  rhsFormOut.close();

  return;
}

void calculateResultFromCache(CustomWeakFormPoisson& wf)
{
  // Utilities.
  int ndof = space->get_num_dofs();
	std::cout << (std::string)"Ndofs: " << ndof << '.' << std::endl;
  Element* e;

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  /// The solution vector.
  double* sln_vector;
	
	if(BASE_VISUALIZATION)
	{
		Views::BaseView<double> b;
		b.show(space, Views::HERMES_EPS_VERYHIGH);
		b.wait_for_close();
	}
	
	Hermes::Mixins::TimeMeasurable time;
	time.tick();

	/// Jacobian.
  SparseMatrix<double>* jacobian = create_matrix<double>();
  jacobian->prealloc(ndof);
  for_all_active_elements(e, mesh)
  {
    AsmList<double> al;
    space->get_element_assembly_list(e, &al);
    for(int i = 0; i < al.get_cnt(); i++)
    {
      if(al.get_dof()[i] < 0)
        continue;
      for(int j = 0; j < al.get_cnt(); j++)
      {
        if(al.get_dof()[j] < 0)
          continue;
        jacobian->pre_add_ij(al.get_dof()[i], al.get_dof()[j]);
        jacobian->pre_add_ij(al.get_dof()[j], al.get_dof()[i]);
      }
    }
  }
  jacobian->alloc();

  /// Residual.
  Vector<double>* residual = create_vector<double>();
  residual->alloc(ndof);
  
  /// Linear solver.
  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);

	int numElements = mesh->get_num_elements();
	int i;

#define CHUNKSIZE 1
	int num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);
#pragma omp parallel private(i) num_threads(num_threads_used)
	{
#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(i = 0; i < numElements; i++)
		{
			Element* e = mesh->get_element(i);
			if(!e->active)
				continue;
			RefMap refmap;
			refmap.set_active_element(e);
			double inv_ref_map_determinant = refmap.get_const_jacobian();
			//double inv_ref_map_determinant = 1 / Hermes::pow(2, P_INIT * 2);
			AsmListEdgeOrientation<double> al;
			space->get_element_assembly_list_with_edge_numbers(e, &al);
			for(int i = 0; i < al.get_cnt(); i++)
			{
				if(al.get_dof()[i] < 0)
					continue;
				for(int j = 0; j < al.get_cnt(); j++)
				{
					if(al.get_dof()[j] >= 0)
					{
						jacobian->add(al.get_dof()[i], al.get_dof()[j], A_values[al.get_idx()[i]][al.get_idx()[j]] * al.get_coef()[i] * al.get_coef()[j]);
					}
					else
					{
						residual->add(al.get_dof()[i], -handleDirichlet(wf, al.get_idx()[i], al.get_idx()[j], e) * al.get_coef()[i] * al.get_coef()[j]);
					}
				}

				residual->add(al.get_dof()[i], rhs_values[al.get_idx()[i]] * al.get_coef()[i] * inv_ref_map_determinant);
			}
		}
	}

	time.tick();
	std::cout << (std::string)"Assembling using cache: " << time.last() << '.' << std::endl;

	try
	{
		FILE* matrixFile = fopen("matrix", "w");
		FILE* rhsFile = fopen("rhs", "w");
		jacobian->dump(matrixFile, "A");
		residual->dump(rhsFile, "b");
		matrix_solver->solve();
	}
	catch(std::exception& e)
	{
		std::cout << e.what();
	}

	time.tick();
	std::cout << (std::string)"Solving: " << time.last() << '.' << std::endl;

  sln_vector = matrix_solver->get_sln_vector();

  // Translate the solution vector into the previously initialized Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, &sln, true);

  // VTK output.
  if(VTK_VISUALIZATION)
  {
    // Output solution in VTK format.
    Hermes::Hermes2D::Views::Linearizer lin;
    bool mode_3D = false;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D, 1, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

    // Output mesh and element orders in VTK format.
    Hermes::Hermes2D::Views::Orderizer ord;
    ord.save_mesh_vtk(space, "mesh.vtk");
    ord.save_orders_vtk(space, "ord.vtk");
  }

  // Visualize the solution.
  Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

  if(HERMES_VISUALIZATION)
  {
    viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_VERYHIGH);
    viewS.wait_for_close();
  }

  return;
}

void calculateResultAssembling(CustomWeakFormPoisson& wf)
{
  // Utilities.
  int ndof = space->get_num_dofs();

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  /// The solution vector.
  double* sln_vector;

	Hermes::Mixins::TimeMeasurable time;
	time.tick();

  SparseMatrix<double>* jacobian = create_matrix<double>();
  Vector<double>* residual = create_vector<double>();
	DiscreteProblemLinear<double> dp(&wf, space);
	dp.assemble(jacobian, residual);

	time.tick();
	std::cout << (std::string)"Ndofs: " << ndof << '.' << std::endl;
	std::cout << (std::string)"Assembling WITHOUT cache: " << time.last() << '.' << std::endl;

  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);
	FILE* matrixFile = fopen("matrix", "w");
	FILE* rhsFile = fopen("rhs", "w");
	jacobian->dump(matrixFile, "A");
	residual->dump(rhsFile, "b");
  fclose(matrixFile);
  fclose(rhsFile);
	matrix_solver->solve();
  sln_vector = matrix_solver->get_sln_vector();

  // Translate the solution vector into the previously initialized Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, &sln, true);

  // VTK output.
  if(VTK_VISUALIZATION)
  {
    // Output solution in VTK format.
    Hermes::Hermes2D::Views::Linearizer lin;
    bool mode_3D = false;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D, 1, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

    // Output mesh and element orders in VTK format.
    Hermes::Hermes2D::Views::Orderizer ord;
    ord.save_mesh_vtk(space, "mesh.vtk");
    ord.save_orders_vtk(space, "ord.vtk");
  }

  // Visualize the solution.
  Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

  if(HERMES_VISUALIZATION)
  {
    viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_VERYHIGH);
    viewS.wait_for_close();
  }

  return;
}

double handleDirichlet(CustomWeakFormPoisson& wf, int shape_indexBasis, int shape_indexDirichlet, Element* e)
{
  RefMap refmap;
  refmap.set_active_element(e);
  Geom<double>* geometry;
  double* jacobian_x_weights;
  int n_quadrature_points = DiscreteProblem<double>::init_geometry_points(&refmap, P_INIT, geometry, jacobian_x_weights);
  PrecalcShapeset pssBasis(space->get_shapeset());
  pssBasis.set_active_element(e);
  PrecalcShapeset pssDirichlet(space->get_shapeset());
  pssDirichlet.set_active_element(e);
  pssBasis.set_active_shape(shape_indexBasis);
  pssDirichlet.set_active_shape(shape_indexDirichlet);
  Func<double>* fnBasis = init_fn(&pssBasis, &refmap, P_INIT);
  Func<double>* fnDirichlet = init_fn(&pssDirichlet, &refmap, P_INIT);
  
  return wf.get_mfvol()[0]->value(n_quadrature_points, jacobian_x_weights, NULL, fnBasis, fnDirichlet, geometry, NULL);
}