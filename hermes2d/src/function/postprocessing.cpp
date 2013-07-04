#include "postprocessing.h"
#include "../space/space.h"
#include "../function/solution.h"
#include <limits>

namespace Hermes
{
  namespace Hermes2D
  {
    namespace PostProcessing
    {
      template<typename Scalar>
      Limiter<Scalar>::Limiter(SpaceSharedPtr<Scalar> space, Scalar* solution_vector) : component_count(1)
      {
        spaces.push_back(space);
        this->init(solution_vector);
      }

      template<typename Scalar>
      Limiter<Scalar>::Limiter(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Scalar* solution_vector) : spaces(spaces), component_count(spaces.size())
      {
        this->init(solution_vector);
      }

      template<typename Scalar>
      void Limiter<Scalar>::init(Scalar* solution_vector_)
      {
        try
        {
          int ndof = Space<Scalar>::get_num_dofs(this->spaces);
          Scalar value = solution_vector_[ndof - 1];

          this->solution_vector = new Scalar[Space<Scalar>::get_num_dofs(this->spaces)];
          memcpy(this->solution_vector, solution_vector_, sizeof(Scalar) * Space<Scalar>::get_num_dofs(this->spaces));
        }
        catch (...)
        {
          throw Exceptions::Exception("Wrong combination of space(s) and solution_vector passed to Limiter().");
        }
      }

      template<typename Scalar>
      MeshFunctionSharedPtr<Scalar> Limiter<Scalar>::get_solution()
      {
        // A check.
        warn_if(this->component_count > 1, "One solution asked from a Limiter, but multiple solutions exist for limiting.");

        return this->get_solutions().back();
      }

      template<typename Scalar>
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > Limiter<Scalar>::get_solutions()
      {
        // Processing.
        this->process();

        if(this->limited_solutions.empty() || this->limited_solutions.size() != this->component_count)
          throw Exceptions::Exception("Limiter failed for unknown reason.");
        else
          return this->limited_solutions;
      }

      template<typename Scalar>
      Hermes::vector<int> Limiter<Scalar>::get_changed_element_ids() const
      {
        return this->changed_element_ids;
      }

      template<typename Scalar>
      bool Limiter<Scalar>::isOkay() const
      {
        bool okay = true;

        return true;
      }

      VertexBasedLimiter::VertexBasedLimiter(SpaceSharedPtr<double> space, double* solution_vector, int maximum_polynomial_order)
        : Limiter<double>(space, solution_vector)
      {
        this->init(maximum_polynomial_order);
      }

      VertexBasedLimiter::VertexBasedLimiter(Hermes::vector<SpaceSharedPtr<double> > spaces, double* solution_vector, int maximum_polynomial_order)
        : Limiter<double>(spaces, solution_vector)
      {
        this->init(maximum_polynomial_order);
      }

      void VertexBasedLimiter::init(int maximum_polynomial_order)
      {
        // Checking that this is the Taylor shapeset.
        for(int i = 0; i < this->component_count; i++)
        {
          if(this->spaces[i]->get_shapeset()->get_id() != 31)
            throw Exceptions::Exception("VertexBasedLimiter designed for L2ShapesetTaylor. Ignore this exception for unforeseen problems.");
        }

        vertex_min_values = NULL;
        vertex_max_values = NULL;

        this->mixed_derivatives_count = (maximum_polynomial_order)*(maximum_polynomial_order + 1) / 2;
      }

      VertexBasedLimiter::~VertexBasedLimiter()
      {
        deallocate_vertex_values();
      }

      void VertexBasedLimiter::process()
      {
        // Start by creating temporary solutions.  
        for(int component = 0; component < this->component_count; component++)
          this->limited_solutions.push_back(MeshFunctionSharedPtr<double>(new Solution<double>(this->spaces[component]->get_mesh())));
        Solution<double>::vector_to_solutions(this->solution_vector, this->spaces, this->limited_solutions);

        // Prepare the vertex values.
        prepare_min_max_vertex_values();

        // Use those to incorporate the correction factor.
        Element* e;
        for(int component = 0; component < this->component_count; component++)
        {
          MeshSharedPtr mesh = this->spaces[component]->get_mesh();

          for_all_active_elements(e, mesh)
          {
            if(this->impose_quadratic_correction_factor(e, component))
              this->impose_linear_correction_factor(e, component);
          }
        }

        // Create the final solutions.
        Solution<double>::vector_to_solutions(this->solution_vector, this->spaces, this->limited_solutions);
      }

      void VertexBasedLimiter::impose_linear_correction_factor(Element* e, int component)
      {
        double correction_factor = std::numeric_limits<double>::infinity();

        double centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, 0);

        AsmList<double> al;
        this->spaces[component]->get_element_assembly_list(e, &al);
        for(int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
        {
          Node* vertex = e->vn[i_vertex];
          double x, y;
          RefMap::untransform(e, vertex->x, vertex->y, x, y);

          double vertex_value = 0.;
          for(int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
            vertex_value += this->solution_vector[al.dof[i_basis_fn]] * this->spaces[component]->get_shapeset()->get_fn_value(al.idx[i_basis_fn], x, y, 0, e->get_mode());

          double fraction;
          if(vertex_value > centroid_value_multiplied)
            fraction = std::min(1., (this->vertex_max_values[component][vertex->id][0] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
          else
            if(vertex_value == centroid_value_multiplied)
              fraction = 1.;
            else
              fraction = std::min(1., (this->vertex_min_values[component][vertex->id][0] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));

          correction_factor = std::min(correction_factor, fraction);
        }
        for(int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
        {
          int order = this->spaces[component]->get_shapeset()->get_order(al.idx[i_basis_fn], e->get_mode());
          if(H2D_GET_H_ORDER(order) == 1 || H2D_GET_V_ORDER(order) == 1)
            this->solution_vector[al.dof[i_basis_fn]] *= correction_factor;
        }
      }

      bool VertexBasedLimiter::impose_quadratic_correction_factor(Element* e, int component)
      {
        double correction_factor = std::numeric_limits<double>::infinity();
        AsmList<double> al;
        this->spaces[component]->get_element_assembly_list(e, &al);

        for(int i_derivative = 1; i_derivative <= 2; i_derivative++)
        {
          double centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, i_derivative);

          for(int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
          {
            Node* vertex = e->vn[i_vertex];
            double x, y;
            RefMap::untransform(e, vertex->x, vertex->y, x, y);

            double vertex_value = 0.;
            for(int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
              vertex_value += this->solution_vector[al.dof[i_basis_fn]] * this->spaces[component]->get_shapeset()->get_value(i_derivative, al.idx[i_basis_fn], x, y, 0, e->get_mode());

            double fraction;
            if(vertex_value > centroid_value_multiplied)
              fraction = std::min(1., (this->vertex_max_values[component][vertex->id][i_derivative] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
            else
              if(vertex_value == centroid_value_multiplied)
                fraction = 1.;
              else
                fraction = std::min(1., (this->vertex_min_values[component][vertex->id][i_derivative] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));

            correction_factor = std::min(correction_factor, fraction);
          }
        }
        for(int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
        {
          int order = this->spaces[component]->get_shapeset()->get_order(al.idx[i_basis_fn], e->get_mode());
          if(H2D_GET_H_ORDER(order) == 2 || H2D_GET_V_ORDER(order) == 2)
            this->solution_vector[al.dof[i_basis_fn]] *= correction_factor;
        }
        if(correction_factor < (1 - 1e-3))
          return true;
        else
          return false;
      }

      void VertexBasedLimiter::prepare_min_max_vertex_values()
      {
        // Reallocate.
        deallocate_vertex_values();
        allocate_vertex_values();

        // Calculate min/max vertex values.
        Element* e;
        for(int component = 0; component < this->component_count; component++)
        {
          MeshSharedPtr mesh = this->spaces[component]->get_mesh();

          for_all_active_elements(e, mesh)
          {
            for(int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
            {
              Node* vertex = e->vn[i_vertex];
              double* element_centroid_value_multiplied = new double[this->mixed_derivatives_count];
              for(int i_derivative = 0; i_derivative < this->mixed_derivatives_count; i_derivative++)
              {
                double element_centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, i_derivative);
                this->vertex_min_values[component][vertex->id][i_derivative] = std::min(this->vertex_min_values[component][vertex->id][i_derivative], element_centroid_value_multiplied);
                this->vertex_max_values[component][vertex->id][i_derivative] = std::max(this->vertex_max_values[component][vertex->id][i_derivative], element_centroid_value_multiplied);
              }
            }
          }
        }
      }

      double VertexBasedLimiter::get_centroid_value_multiplied(Element* e, int component, int mixed_derivative_index)
      {
        if(mixed_derivative_index > 5)
        {
          throw Exceptions::MethodNotImplementedException("VertexBasedLimiter::get_centroid_value_multiplied only works for first and second derivatives.");
          return 0.;
        }

        Solution<double>* sln = dynamic_cast<Solution<double>*>(this->limited_solutions[component].get());
        double result;
        if(e->get_mode() == HERMES_MODE_TRIANGLE)
          result = sln->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, mixed_derivative_index);
        else
          result = sln->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, mixed_derivative_index);

        switch(mixed_derivative_index)
        {
        case 1:
            result *= ELEMENT_DELTA_X;
          break;
        case 2:
            result *= ELEMENT_DELTA_Y;
          break;
        case 3:
            result *= ELEMENT_DELTA_X;
            result *= ELEMENT_DELTA_X;
          break;
        case 4:
            result *= ELEMENT_DELTA_Y;
            result *= ELEMENT_DELTA_Y;
          break;
        case 5:
            result *= ELEMENT_DELTA_Y;
            result *= ELEMENT_DELTA_X;
          break;
        }
      }

      void VertexBasedLimiter::allocate_vertex_values()
      {
        this->vertex_min_values = new double**[this->component_count];
        for(int i = 0; i < this->component_count; i++)
        {
          this->vertex_min_values[i] = new double*[this->spaces[i]->get_mesh()->get_max_node_id()];

          for(int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
          {
            this->vertex_min_values[i][j] = new double[this->spaces[i]->get_mesh()->get_max_node_id()];
            for(int k = 0; k < this->mixed_derivatives_count; k++)
              this->vertex_min_values[i][j][k] = std::numeric_limits<double>::infinity();
          }
        }

        this->vertex_max_values = new double**[this->component_count];
        for(int i = 0; i < this->component_count; i++)
        {
          this->vertex_max_values[i] = new double*[this->spaces[i]->get_mesh()->get_max_node_id()];

          for(int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
          {
            this->vertex_max_values[i][j] = new double[this->spaces[i]->get_mesh()->get_max_node_id()];
            for(int k = 0; k < this->mixed_derivatives_count; k++)
              this->vertex_max_values[i][j][k] = -std::numeric_limits<double>::infinity();
          }
        }
      }

      void VertexBasedLimiter::deallocate_vertex_values()
      {
        if(this->vertex_min_values)
        {
          for(int i = 0; i < this->component_count; i++)
          {
            for(int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
            {
              delete [] this->vertex_min_values[i][j];
              delete [] this->vertex_max_values[i][j];
            }

            delete [] this->vertex_min_values[i];
            delete [] this->vertex_max_values[i];
          }

          delete [] this->vertex_min_values;
          delete [] this->vertex_max_values;
        }
      }

      template class HERMES_API Limiter<double>;
      template class HERMES_API Limiter<std::complex<double> >;
    }
  }
}