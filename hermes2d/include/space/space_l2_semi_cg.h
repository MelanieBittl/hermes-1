#ifndef __H2D_SPACE_L2_SEMI_CG
#define __H2D_SPACE_L2_SEMI_CG

#include "space.h"
#include "../shapeset/shapeset_l2_semi.h"
#include "../shapeset/shapeset_bb.h"
namespace Hermes
{
namespace Hermes2D
{
/// @ingroup spaces

///L2_SEMI_CG_SPACE represents a space of continuous linear functions enriched with discontinuous basisfunctions for higher orders..
template<typename Scalar>
class HERMES_API L2_SEMI_CG_Space : public Space<Scalar>
{
public:
L2_SEMI_CG_Space(bool serendipity= false);
L2_SEMI_CG_Space(MeshSharedPtr mesh, EssentialBCs<Scalar>* boundary_conditions, int p_init = 1,bool serendipity= false,
Shapeset* shapeset = NULL);

L2_SEMI_CG_Space(MeshSharedPtr mesh, int p_init = 1,bool serendipity = false,Shapeset* shapeset = NULL);

virtual ~L2_SEMI_CG_Space();

virtual void set_shapeset(Shapeset* shapeset);

      /// Obtains an assembly list for the given element.
      virtual void get_element_assembly_list(Element* e, AsmList<Scalar>* al) const;

/// Removes the degree of freedom from a vertex node with the given id (i.e., its number
/// in the mesh file) and makes it part of the Dirichlet lift with the given value.
/// This is a special-purpose function which normally should not be needed.
/// It is intended for fixing the solution of a system which would otherwise be singular
/// and for some reason a standard Dirichlet condition (with non-zero measure on the
/// boundary) is not suitable.
void fix_vertex(int id, Scalar value = 0.0);

virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc);

/// Copy from Space instance 'space'
virtual void copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh);

protected:

virtual int get_edge_order(Element* e, int edge) const {
return H2D_MAKE_EDGE_ORDER(e->get_mode(), edge, this->edata[e->id].order);
}

bool serendipity; 
virtual SpaceType get_type() const { return HERMES_L2_SEMI_SPACE; }

/// Common code for the constructors.
void init(Shapeset* shapeset, int p_init);

virtual void assign_vertex_dofs();	
virtual void assign_edge_dofs(){};	
virtual void assign_bubble_dofs();

      /// Obtains an edge assembly list 
      virtual void get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al) const;

virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const;
virtual void get_boundary_assembly_list_internal(Element* e, int ie, AsmList<Scalar>* al) const;
virtual void get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const;


struct EdgeInfo
{
Node* node;
int part;
int ori;
double lo, hi;
};


      virtual void update_edge_bc(Element* e, SurfPos* surf_pos);

inline void output_component(typename Space<Scalar>::BaseComponent*& current, typename Space<Scalar>::BaseComponent*& last, typename Space<Scalar>::BaseComponent* min);

typename Space<Scalar>::BaseComponent* merge_baselists(typename Space<Scalar>::BaseComponent* l1, int n1, typename Space<Scalar>::BaseComponent* l2, int n2, int& ncomponents);

void update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3);

virtual void update_constraints();

friend class Space<Scalar>;

};
}
}









#endif
