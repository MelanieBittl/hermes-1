// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_SOLUTION_H
#define __H2D_SOLUTION_H

#include "../function/function.h"
#include "../space/space.h"
#include "../mesh/refmap.h"
#include "../../../hermes_common/matrix.h"

class PrecalcShapeset;
class Ord;

/// \brief Represents a function defined on a mesh.
///
/// MeshFunction is a base class for all classes representing an arbitrary function
/// superimposed on a mesh (ie., domain). These include the Solution, ExactSolution
/// and Filter classes, which define the concrete behavior and the way the function
/// is (pre)calculated. Any such function can later be visualized.
///
/// (This is an abstract class and cannot be instantiated.)
///
template<typename Scalar>
class HERMES_API MeshFunction : public Function<Scalar>
{
public:

  MeshFunction();
  MeshFunction(Mesh *mesh);
  virtual ~MeshFunction();

  virtual void init() {};
  virtual void reinit() { this->free(); init();};

  virtual void set_quad_2d(Quad2D* quad_2d);
  virtual void set_active_element(Element* e);

  virtual int get_edge_fn_order(int edge) { return get_edge_fn_order(edge); }

  Mesh*   get_mesh() const { return mesh; }
  RefMap* get_refmap() { this->update_refmap(); return refmap; }

  virtual Scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0) = 0;

  /// Virtual function handling overflows. Has to be virtual, because
  /// the necessary iterators in the templated class do not work with GCC.
  virtual void handle_overflow_idx();

  /// See Transformable::push_transform.
  virtual void push_transform(int son);

  virtual void pop_transform();
  
protected:

  int mode;
  Mesh* mesh;
  RefMap* refmap;

public:

  /// For internal use only.
  void force_transform(MeshFunction<Scalar>* mf)
    { Function<Scalar>::force_transform(mf->get_transform(), mf->get_ctm()); }
  void update_refmap()
    { refmap->force_transform(this->sub_idx, this->ctm); }
  void force_transform(uint64_t sub_idx, Trf* ctm)
  {
    this->sub_idx = sub_idx;
    this->ctm = ctm;
  }
};

/// \brief Represents the solution of a PDE.
///
/// The Solution class represents the solution of a PDE. Given a space and a solution vector,
/// it calculates the appropriate linear combination of basis functions at the specified
/// element and integration points.
///
/// TODO: write how to obtain solution values, maybe include inherited methods from Function as comments.
///
template<typename Scalar>
class HERMES_API Solution : public MeshFunction<Scalar>
{
public:

  void init();
  Solution();
  Solution(Mesh *mesh);
  Solution(Mesh *mesh, Scalar init_const);
  Solution(Mesh *mesh, Scalar init_const_0, Scalar init_const_1);
  Solution (Space<Scalar>* s, Vector<Scalar>* coeff_vec);
  Solution (Space<Scalar>* s, Scalar* coeff_vec);
  virtual ~Solution();
  virtual void free();

  void assign(Solution<Scalar>* sln);
  Solution& operator = (Solution& sln) { assign(&sln); return *this; }
  void copy(const Solution<Scalar>* sln);

  int* get_element_orders() { return this->elem_orders;}

  void set_const(Mesh* mesh, Scalar c);
  void set_const(Mesh* mesh, Scalar c0, Scalar c1); // two-component (Hcurl) const

  void set_zero(Mesh* mesh);
  void set_zero_2(Mesh* mesh); // two-component (Hcurl) zero

  virtual int get_edge_fn_order(int edge) { return MeshFunction<Scalar>::get_edge_fn_order(edge); }
  int get_edge_fn_order(int edge, Space<Scalar>* space, Element* e = NULL);

  /// Sets solution equal to Dirichlet lift only, solution vector = 0
  void set_dirichlet_lift(Space<Scalar>* space, PrecalcShapeset* pss = NULL);

  /// Enables or disables transformation of the solution derivatives (H1 case)
  /// or values (vector (Hcurl) case). This means H2D_FN_DX_0 and H2D_FN_DY_0 or
  /// H2D_FN_VAL_0 and H2D_FN_VAL_1 will or will not be returned premultiplied by the reference
  /// mapping matrix. The default is enabled (true).
  void enable_transform(bool enable = true);

  /// Saves the complete solution (i.e., including the internal copy of the mesh and
  /// element orders) to a binary file. On Linux, if `compress` is true, the file is
  /// compressed with gzip and a ".gz" suffix added to the file name.
  void save(const char* filename, bool compress = true);

  /// Loads the solution from a file previously created by Solution::save(). This completely
  /// restores the solution in the memory. The file name has to include the ".gz" suffix,
  /// in which case the file is piped through gzip to decompress the data (Linux only).
  void load(const char* filename);

  /// Returns solution value or derivatives at element e, in its reference domain point (xi1, xi2).
  /// 'item' controls the returned value: 0 = value, 1 = dx, 2 = dy, 3 = dxx, 4 = dyy, 5 = dxy.
  /// NOTE: This function should be used for postprocessing only, it is not effective
  /// enough for calculations.
  Scalar get_ref_value(Element* e, double xi1, double xi2, int component = 0, int item = 0);

  /// Returns solution value or derivatives (correctly transformed) at element e, in its reference
  /// domain point (xi1, xi2). 'item' controls the returned value: 0 = value, 1 = dx, 2 = dy,
  /// 3 = dxx, 4 = dyy, 5 = dxy.
  /// NOTE: This function should be used for postprocessing only, it is not effective
  /// enough for calculations.
  Scalar get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b);

  /// Returns solution value or derivatives at the physical domain point (x, y).
  /// 'item' controls the returned value: H2D_FN_VAL_0, H2D_FN_VAL_1, H2D_FN_DX_0, H2D_FN_DX_1, H2D_FN_DY_0,....
  /// NOTE: This function should be used for postprocessing only, it is not effective
  /// enough for calculations. Since it searches for an element sequentinally, it is extremelly
  /// slow. Prefer Solution::get_ref_value if possible.
  virtual Scalar get_pt_value(double x, double y, int item = H2D_FN_VAL_0);

  /// Returns the number of degrees of freedom of the solution.
  /// Returns -1 for exact or constant solutions.
  int get_num_dofs() const { return num_dofs; };

  /// Multiplies the function represented by this class by the given coefficient.
  void multiply(Scalar coef);

  /// Returns solution type.
  ESolutionType get_type() const { return sln_type; };

  /// Returns space type.
  ESpaceType get_space_type() const { return space_type; };

public:
  /// Internal.
  virtual void set_active_element(Element* e);

  /// Passes solution components calculated from solution vector as Solutions.
  static void vector_to_solutions(Scalar* solution_vector, Hermes::vector<Space<Scalar>*> spaces,
                                  Hermes::vector<Solution<Scalar>*> solutions,
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(Scalar* solution_vector, Space<Scalar>* space, Solution<Scalar>* solution,
                                 bool add_dir_lift = true);
  static void vector_to_solutions(Vector<Scalar>* vec, Hermes::vector<Space<Scalar>*> spaces,
                                  Hermes::vector<Solution<Scalar>*> solutions,
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(Vector<Scalar>* vec, Space<Scalar>* space, Solution<Scalar>* solution,
                                 bool add_dir_lift = true);
  static void vector_to_solutions(Scalar* solution_vector, Hermes::vector<Space<Scalar>*> spaces,
                                  Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<PrecalcShapeset *> pss,
                                  Hermes::vector<bool> add_dir_lift = Hermes::vector<bool>());
  static void vector_to_solution(Scalar* solution_vector, Space<Scalar>* space, Solution<Scalar>* solution,
                                 PrecalcShapeset* pss, bool add_dir_lift = true);

  bool own_mesh;
protected:

  /// Converts a coefficient vector into a Solution.
  virtual void set_coeff_vector(Space<Scalar>* space, Vector<Scalar>* vec, bool add_dir_lift);
  virtual void set_coeff_vector(Space<Scalar>* space, PrecalcShapeset* pss, Scalar* coeffs, bool add_dir_lift);
  virtual void set_coeff_vector(Space<Scalar>* space, Scalar* coeffs, bool add_dir_lift);

  ESolutionType sln_type;

  bool transform;

  /// Precalculated tables for last four used elements.
  /// There is a 2-layer structure of the precalculated tables.
  /// The first (the lowest) one is the layer where mapping of integral orders to
  /// Function::Node takes place. See function.h for details.
  /// The second one is the layer with mapping of sub-element transformation to
  /// a table from the lowest layer.
  /// The highest layer (in contrast to the PrecalcShapeset class) is represented
  /// here only by this array.
  std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>* tables[4][4];

  Element* elems[4][4];
  int cur_elem, oldest[4];

  Scalar* mono_coefs;  ///< monomial coefficient array
  int* elem_coefs[2];  ///< array of pointers into mono_coefs
  int* elem_orders;    ///< stored element orders
  int num_coefs, num_elems;
  int num_dofs;

  ESpaceType space_type;
  void transform_values(int order, struct Function<Scalar>::Node* node, int newmask, int oldmask, int np);

  Scalar   cnst[2];
  Scalar   exact_mult;

  virtual void precalculate(int order, int mask);

  Scalar* dxdy_coefs[2][6];
  Scalar* dxdy_buffer;

  double** calc_mono_matrix(int o, int*& perm);
  void init_dxdy_buffer();
  void free_tables();

  Element* e_last; ///< last visited element when getting solution values at specific points

};


/// \brief Represents an exact solution of a PDE.
///
/// ExactSolution represents an arbitrary user-specified function defined on a domain (mesh),
/// typically an exact solution to a PDE. This can be used to compare an approximate solution
/// with an exact solution (see DiffFilter).
template<typename Scalar>
class ExactSolution : public Solution<Scalar>
{
public:
  ExactSolution(Mesh* mesh);

  ~ExactSolution();

  // Dimension of result - either 1 or 2.
  virtual unsigned int get_dimension() const = 0;
};

/// These classes are abstract (pure virtual destructor).
/// The user is supposed to subclass them (see e.g. NIST benchmarks).
template<typename Scalar>
class HERMES_API ExactSolutionScalar : public ExactSolution<Scalar>
{
public:
  ExactSolutionScalar(Mesh* mesh);

  ~ExactSolutionScalar() = 0;

  // For Scalar-valued solutions this returns 1.
  virtual unsigned int get_dimension() const;

  // Function returning the value.
  virtual Scalar value (double x, double y) const = 0;

  // Function returning the derivatives.
  virtual void derivatives (double x, double y, Scalar& dx, Scalar& dy) const = 0;

  // Function returning the value and derivatives.
  Scalar exact_function (double x, double y, Scalar& dx, Scalar& dy) const {
    derivatives (x, y, dx, dy);
    return value (x, y);
  };

  // Function returning the integration order that 
  // should be used when integrating the function.
  virtual Ord ord(Ord x, Ord y) const = 0;
};

template<typename Scalar>
class HERMES_API ExactSolutionVector : public ExactSolution<Scalar>
{
public:
  ExactSolutionVector(Mesh* mesh);

  ~ExactSolutionVector() = 0;

  // For vector-valued solutions this returns 2.
  virtual unsigned int get_dimension() const;

  // Function returning the value.
  virtual Scalar2<Scalar> value (double x, double y) const = 0;

  // Function returning the derivatives.
  virtual void derivatives (double x, double y, Scalar2<Scalar>& dx, Scalar2<Scalar>& dy) const = 0;

  // Function returning the value and derivatives.
  virtual Scalar2<Scalar> exact_function(double x, double y, Scalar2<Scalar>& dx, Scalar2<Scalar>& dy) const {
    derivatives (x, y, dx, dy);
    return value (x, y);
  };

  // Function returning the integration order that 
  // should be used when integrating the function.
  virtual Ord ord(Ord x, Ord y) const = 0;
};
#endif
