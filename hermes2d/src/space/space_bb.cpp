#include "space_bb.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    SpaceBB<Scalar>::SpaceBB() : Space<Scalar>()
    {
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::init(Shapeset* shapeset, int p_init)
    {
      if(shapeset == NULL)
      {
        this->shapeset = new ShapesetBB(p_init);
				//this->shapeset = new H1Shapeset;
        this->own_shapeset = true;
      }

     // this->precalculate_projection_matrix(2, this->proj_mat, this->chol_p);

      // set uniform poly order in elements
      if(p_init < 1) 
        throw Hermes::Exceptions::Exception("P_INIT must be >=  1 in an H1 space.");

      else this->set_uniform_order_internal(p_init, HERMES_ANY_INT);

      // enumerate basis functions
      this->assign_dofs();
    }


    template<typename Scalar>
    SpaceBB<Scalar>::SpaceBB(MeshSharedPtr mesh, int p_init, Shapeset* shapeset)
      : Space<Scalar>(mesh, shapeset, NULL)
    {
      init(shapeset, p_init);
    }

    template<typename Scalar>
    SpaceBB<Scalar>::~SpaceBB()
    {
      if(this->own_shapeset)
        delete this->shapeset;
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh)
    {
      Space<Scalar>::copy(space, new_mesh);

     // this->precalculate_projection_matrix(2, this->proj_mat, this->chol_p);

      this->assign_dofs();
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::set_shapeset(Shapeset *shapeset)
    {
      if(shapeset->get_id() ==5)
      {
        this->shapeset = shapeset;
        this->own_shapeset = false;
      }
      else
        throw Hermes::Exceptions::Exception("Wrong shapeset type in SpaceBB<Scalar>::set_shapeset()");
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::assign_vertex_dofs()
    {

      // Vertex dofs.
      Element* e;
      this->vertex_functions_count = 0;
      for_all_active_elements(e, this->mesh)
      {
        int order = this->get_element_order(e->id);
        if(order > 0)
        {
          for (unsigned int i = 0; i < e->get_nvert(); i++)
          {
            Node* vn = e->vn[i];
            typename Space<Scalar>::NodeData* nd = this->ndata + vn->id;
            if(!vn->is_constrained_vertex() && nd->dof == this->H2D_UNASSIGNED_DOF)
            {
               nd->dof = this->next_dof;
                this->next_dof += this->stride;
                this->vertex_functions_count++;
              
              nd->n = 1;
            }
          }
        }
      }
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::assign_edge_dofs()
    {
      // Edge dofs.
      Element* e;
      this->edge_functions_count = 0;
      for_all_active_elements(e, this->mesh)
      {
        int order = this->get_element_order(e->id);
        if(order > 0)
        {
          for (unsigned int i = 0; i < e->get_nvert(); i++)
          {
            Node* vn = e->vn[i];
            typename Space<Scalar>::NodeData* nd = this->ndata + vn->id;
            Node* en = e->en[i];
            nd = this->ndata + en->id;
            if(nd->dof == this->H2D_UNASSIGNED_DOF)
            {
              // If the edge node is not constrained, assign it dofs.
              if(en->ref > 1 || en->bnd || this->mesh->peek_vertex_node(en->p1, en->p2) != NULL)
              {
                int ndofs = this->get_edge_order_internal(en) - 1;
                nd->n = ndofs;
                  nd->dof = this->next_dof;
                  this->next_dof += ndofs * this->stride;
                  this->edge_functions_count += ndofs;
                
              }
              else // Constrained edge node.
                nd->n = -1;
            }
          }
        }
      }

    }
   template<typename Scalar>
    void SpaceBB<Scalar>::assign_bubble_dofs()
    {
      // Bubble dofs.
      Element* e;
      this->bubble_functions_count = 0;
      for_all_active_elements(e, this->mesh)
      {
        int order = this->get_element_order(e->id);
        if(order > 0)
        {
          typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
          ed->bdof = this->next_dof;
          ed->n = this->shapeset->get_num_bubbles(ed->order, e->get_mode());
          this->next_dof += ed->n * this->stride;
          this->bubble_functions_count += ed->n;
        }
      }
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const
    {
      Node* vn = e->vn[iv];
      typename Space<Scalar>::NodeData* nd = &this->ndata[vn->id];
      int index = this->shapeset->get_vertex_index(iv, e->get_mode());
      if(this->get_element_order(e->id) == 0) return;

      if(!vn->is_constrained_vertex()) // unconstrained
      {
        al->add_triplet(index, nd->dof, (nd->dof >= 0) ? 1.0 : *(nd->vertex_bc_coef));
      }
      else // constrained
      {
        //debug_log("! B cause of the triplet\n");
        for (int j = 0; j < nd->ncomponents; j++)
          if(nd->baselist[j].coef != (Scalar) 0)
          {
            al->add_triplet(index, nd->baselist[j].dof, nd->baselist[j].coef);
          }
      }
    }

    template<typename Scalar>
    void SpaceBB<Scalar>::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const
    {
      Node* en = e->en[surf_num];
      typename Space<Scalar>::NodeData* nd = &this->ndata[en->id];
      if(this->get_element_order(e->id) == 0)
        return;

      if(nd->n >= 0) // unconstrained
      {
        if(nd->dof >= 0)
        {
          int ori = (e->vn[surf_num]->id < e->vn[e->next_vert(surf_num)]->id) ? 0 : 1;
          for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += this->stride)
            al->add_triplet(this->shapeset->get_edge_index(surf_num, ori, j + 2, e->get_mode()), dof, 1.0);
        }
        else
        {
          for (int j = 0; j < nd->n; j++)
          {
            al->add_triplet(this->shapeset->get_edge_index(surf_num, 0, j + 2, e->get_mode()), -1, nd->edge_bc_proj[j + 2]);
          }
        }
      }
      else // constrained
      {
        int part = nd->part;
        int ori = part < 0 ? 1 : 0;
        if(part < 0) part ^=  ~0;

        nd = &this->ndata[nd->base->id];
        for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += this->stride)
          al->add_triplet(this->shapeset->get_constrained_edge_index(surf_num, j + 2, ori, part, e->get_mode()), dof, 1.0);
      }
    }



 template<typename Scalar>
    Scalar* SpaceBB<Scalar>::get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc)
    {
      assert(order >= 1);
      Scalar* proj = new Scalar[order + 1];

      if(bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_CONST)
      {
        proj[0] = proj[1] = bc->value_const;
      } // If the BC is not constant.
      else if(bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_FUNCTION)
      {
        surf_pos->t = surf_pos->lo;
        // Find out the (x, y) coordinates for the first endpoint.
        double x, y, n_x, n_y, t_x, t_y;
        Nurbs* nurbs = surf_pos->base->is_curved() ? surf_pos->base->cm->nurbs[surf_pos->surf_num] : NULL;
        CurvMap::nurbs_edge(surf_pos->base, nurbs, surf_pos->surf_num, 2.0*surf_pos->t - 1.0, x, y, n_x, n_y, t_x, t_y);
        // Calculate.
        proj[0] = bc->value(x, y, n_x, n_y, t_x, t_y);
        surf_pos->t = surf_pos->hi;
        // Find out the (x, y) coordinates for the second endpoint.
        CurvMap::nurbs_edge(surf_pos->base, nurbs, surf_pos->surf_num, 2.0*surf_pos->t - 1.0, x, y, n_x, n_y, t_x, t_y);
        // Calculate.
        proj[1] = bc->value(x, y, n_x, n_y, t_x, t_y);
      }

      if(order-- > 1)
      {
        Quad1DStd quad1d;
        Scalar* rhs = proj + 2;
        int mo = quad1d.get_max_order();
        double2* pt = quad1d.get_points(mo);

        // get boundary values at integration points, construct rhs
        for (int i = 0; i < order; i++)
        {
          rhs[i] = 0.0;
          int ii = this->shapeset->get_edge_index(0, 0, i + 2, surf_pos->base->get_mode());
          for (int j = 0; j < quad1d.get_num_points(mo); j++)
          {
            double t = (pt[j][0] + 1) * 0.5, s = 1.0 - t;
            Scalar l = proj[0] * s + proj[1] * t;
            surf_pos->t = surf_pos->lo * s + surf_pos->hi * t;

            if(bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_CONST)
              rhs[i] += pt[j][1] * this->shapeset->get_fn_value(ii, pt[j][0], -1.0, 0, surf_pos->base->get_mode())
              * (bc->value_const - l);
            // If the BC is not constant.
            else if(bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_FUNCTION)
            {
              // Find out the (x, y) coordinate.
              double x, y, n_x, n_y, t_x, t_y;
              Nurbs* nurbs = surf_pos->base->is_curved() ? surf_pos->base->cm->nurbs[surf_pos->surf_num] : NULL;
              CurvMap::nurbs_edge(surf_pos->base, nurbs, surf_pos->surf_num, 2.0*surf_pos->t - 1.0, x, y, n_x, n_y, t_x, t_y);
              // Calculate.
              rhs[i] += pt[j][1] * this->shapeset->get_fn_value(ii, pt[j][0], -1.0, 0, surf_pos->base->get_mode())
                * (bc->value(x, y, n_x, n_y, t_x, t_y) - l);
            }
          }
        }

        // solve the system using a precalculated Cholesky decomposed projection matrix
        cholsl(this->proj_mat, order, this->chol_p, rhs, rhs);
      }

      return proj;
    }

 
    template HERMES_API class SpaceBB<double>;
    template HERMES_API class SpaceBB<std::complex<double> >;
  }
}
