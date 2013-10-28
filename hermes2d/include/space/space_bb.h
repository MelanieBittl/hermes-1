#ifndef __H2D_SPACE_BB_H
#define __H2D_SPACE_BB_H

#include "../shapeset/shapeset_bb.h"
//#include "../shapeset/shapeset_h1_all.h"
#include "space.h"

namespace Hermes
{
	namespace Hermes2D
	{
		template<typename Scalar>
		class HERMES_API SpaceBB : public Space<Scalar>
		{
		public:
			SpaceBB();
			SpaceBB(MeshSharedPtr mesh, int p_init = 1,	Shapeset* shapeset = NULL);
      SpaceBB(MeshSharedPtr mesh, EssentialBCs<Scalar>* boundary_conditions, int p_init = 1,
        Shapeset* shapeset = NULL);

			virtual ~SpaceBB();

			virtual void set_shapeset(Shapeset* shapeset);

			/// Copy from Space instance 'space'
			virtual void copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh);

      virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc);

		protected:

			virtual SpaceType get_type() const { return HERMES_H1_SPACE; }

			/// Common code for the constructors.
			void init(Shapeset* shapeset, int p_init);

			virtual void assign_vertex_dofs();
			virtual void assign_edge_dofs();
			virtual void assign_bubble_dofs();

			virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const;
			virtual void get_boundary_assembly_list_internal(Element* e, int ie, AsmList<Scalar>* al) const;

			struct EdgeInfo
			{
				Node* node;
				int part;
				int ori;
				double lo, hi;
			};


			inline void output_component(typename Space<Scalar>::BaseComponent*& current, typename Space<Scalar>::BaseComponent*& last, typename Space<Scalar>::BaseComponent* min,
				Node*& edge, typename Space<Scalar>::BaseComponent*& edge_dofs, int order);

			typename Space<Scalar>::BaseComponent* merge_baselists(typename Space<Scalar>::BaseComponent* l1, int n1, typename Space<Scalar>::BaseComponent* l2, int n2,
				Node* edge, typename Space<Scalar>::BaseComponent*& edge_dofs, int& ncomponents, int order);

			void update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3);

			virtual void update_constraints();

      virtual void precalculate_projection_matrix(int nv, double**& mat, double*& p);
			friend class Space<Scalar>;
		};
	}
}
#endif
