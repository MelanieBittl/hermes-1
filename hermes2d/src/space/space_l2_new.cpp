 #include "space_l2_new.h"
 
// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D. If not, see <http://www.gnu.org/licenses/>.


namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    L2_NEW_Space<Scalar>::L2_NEW_Space(bool serendipity) : L2Space<Scalar>(), serendipity(serendipity)
    {
			
    }

    template<typename Scalar>
    void L2_NEW_Space<Scalar>::init(Shapeset* shapeset, int p_init)
    {
      if(shapeset == NULL)
      {
        this->shapeset = new H1Shapeset;
        this->own_shapeset = true;
      }
      this->precalculate_projection_matrix(2, this->proj_mat, this->chol_p);

      // set uniform poly order in elements
      if(p_init < 1)
        throw Hermes::Exceptions::Exception("P_INIT must be >= 1 in an L2_NEW_Space space.");

      else this->set_uniform_order_internal(p_init, HERMES_ANY_INT);

      // enumerate basis functions
      this->assign_dofs();
    }


    template<typename Scalar>
    L2_NEW_Space<Scalar>::L2_NEW_Space(MeshSharedPtr mesh, int p_init,bool serendipity, Shapeset* shapeset)
      : L2Space<Scalar>(mesh, p_init , shapeset), serendipity(serendipity)
    {
      init(shapeset, p_init);
    }

    template<typename Scalar>
    L2_NEW_Space<Scalar>::~L2_NEW_Space()
    {
      if(this->own_shapeset)
        delete this->shapeset;
    }

    template<typename Scalar>
    void L2_NEW_Space<Scalar>::copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh)
    {
      Space<Scalar>::copy(space, new_mesh);
      this->precalculate_projection_matrix(2, this->proj_mat, this->chol_p);
      this->assign_dofs();
    }

    template<typename Scalar>
    void L2_NEW_Space<Scalar>::set_shapeset(Shapeset *shapeset)
    {
     if((this->shapeset!=NULL)&& (this->own_shapeset))
     			delete this->shapeset;
			if(shapeset->get_id()<10)
      {
        this->shapeset = shapeset;
        this->own_shapeset = false;
      }
      else
        throw Hermes::Exceptions::Exception("Wrong shapeset type in L2_NEW_Space<Scalar>::set_shapeset()");
    }

    
    
    template<typename Scalar>
    void L2_NEW_Space<Scalar>::assign_bubble_dofs()
    {
      // vertex, Bubble & edge dofs.
      Element* e;
		int ndofs;
      this->bubble_functions_count = 0;
      for_all_active_elements(e, this->mesh)
      {
        int order = this->get_element_order(e->id);
          int ndofs_total = 0;
          typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
          ed->bdof = this->next_dof;
       
		 if(order > 0)
        {
          for (unsigned int i = 0; i < e->get_nvert(); i++)
          {
			ndofs_total++;
          }
        }
        if(order > 1)
        {

          for (unsigned int i = 0; i < e->get_nvert(); i++)
          {
             Node* en = e->en[i];
			if(en->ref >= 1 || en->bnd)
			{
				ndofs = this->get_edge_order_internal(en) - 1;	
				ndofs_total += ndofs;
			}
          }
          
			if(!serendipity)
			{ 		ndofs = this->shapeset->get_num_bubbles(ed->order, e->get_mode()) ;
					ndofs_total += ndofs;
			}

        }
          ed->n = ndofs_total;
          this->next_dof += ed->n * this->stride;
          this->bubble_functions_count += ed->n;
      }
    }

 /*   template<typename Scalar>
    void L2_NEW_Space<Scalar>::get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const
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
        for (int j = 0; j < nd->ncomponents; j++)
          if(nd->baselist[j].coef != (Scalar) 0)
          {
            al->add_triplet(index, nd->baselist[j].dof, nd->baselist[j].coef);
          }
      }
    }*/

    template<typename Scalar>
    void L2_NEW_Space<Scalar>::get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al) const
    {
      this->check();
      al->cnt = 0;
		this->get_element_assembly_list( e, al) ;

    }


    template<typename Scalar>
    void L2_NEW_Space<Scalar>::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const
    { //edge basis functions
        
   /*  Node* en = e->en[surf_num];
      typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
      if(this->get_element_order(e->id) == 0)
        return;
      if(!ed->n) return;
      else
      {
  //int ndofs_total = 0;
  int ndofs_start = 0;int ndofs =0;
  for (unsigned int i = 0; i <= surf_num; i++)
          {
						Node* en = e->en[i];
						if(en->ref >= 1 || en->bnd)
						{
						ndofs = this->get_edge_order_internal(en) - 1;	
						//ndofs_total += ndofs;
						if(i<surf_num) ndofs_start +=ndofs;

						}
            
          }
 
  ndofs_start += ed->bdof;
          int ori = (e->vn[surf_num]->id < e->vn[e->next_vert(surf_num)]->id) ? 0 : 1;
          for (int j = 0, dof = ndofs_start; j < ndofs; j++, dof += this->stride)
            al->add_triplet(this->shapeset->get_edge_index(surf_num, ori, j + 2, e->get_mode()), dof, 1.0);
       }*/
     
    }
    
        template<typename Scalar>
    void L2_NEW_Space<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const
    {

      typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
      if(!ed->n) return;
         int ndofs_total = 0;
			int ndofs;
          int ndofs_start = ed->bdof;
		for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
		{
			 
			  int index = this->shapeset->get_vertex_index(iv, e->get_mode());
			  if(this->get_element_order(e->id) == 0) return;

				al->add_triplet(index, ndofs_start, 1.0 ); 
				ndofs_start += this->stride;
				ndofs_total++;

		  } 


//edge and bubble 
					for (unsigned int i = 0; i < e->get_nvert(); i++)
					{
					Node* en = e->en[i]; int surf_num = i;
					if(en->ref >= 1 || en->bnd )
						{
						ndofs = this->get_edge_order_internal(en) - 1;	
						ndofs_total += ndofs;
						int ori = (e->vn[surf_num]->id < e->vn[e->next_vert(surf_num)]->id) ? 0 : 1;
											 for (int j = 0, dof = ndofs_start; j < ndofs; j++, dof += this->stride)
												 al->add_triplet(this->shapeset->get_edge_index(surf_num, ori, j + 2, e->get_mode()), dof, 1.0);	
						ndofs_start +=ndofs;
						}
					}
				if(!serendipity)
    { int* indices = this->shapeset->get_bubble_indices(ed->order, e->get_mode());
      for (int i = 0, dof = ndofs_start; i < (ed->n - ndofs_total); i++, dof += this->stride)
        al->add_triplet(*indices++, dof, 1.0);
		}
      
    }
    
        template<typename Scalar>
    void L2_NEW_Space<Scalar>::get_element_assembly_list(Element* e, AsmList<Scalar>* al) const
    {
          // some checks
    if(e->id >= this->esize || this->edata[e->id].order < 0)
       throw Hermes::Exceptions::Exception("Uninitialized element order in get_element_assembly_list(id = #%d).", e->id);
    if(!this->is_up_to_date())
        throw Hermes::Exceptions::Exception("The space in get_element_assembly_list() is out of date. You need to update it with assign_dofs()"" any time the mesh changes.");

      // add vertex, edge and bubble functions to the assembly list
      al->cnt = 0;
			get_bubble_assembly_list(e, al);
           
    }
    
   

    template<typename Scalar>
    Scalar* L2_NEW_Space<Scalar>::get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc)
    {
            throw Hermes::Exceptions::Exception("Method get_bc_projection() called from an L2Space.");
      return NULL;
    }


    template HERMES_API class L2_NEW_Space<double>;
    template HERMES_API class L2_NEW_Space<std::complex<double> >;
  }
}
