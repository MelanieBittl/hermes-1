#ifndef __H2D_SPACE_L2_NEW
#define __H2D_SPACE_L2_NEW

#include "space_l2.h"
#include "../shapeset/shapeset_h1_all.h"

namespace Hermes
{
namespace Hermes2D
{
/// @ingroup spaces

///L2_SEMI_CG_SPACE represents a space of continuous linear functions enriched with discontinuous basisfunctions for higher orders..
template<typename Scalar>
class HERMES_API L2_NEW_Space : public L2Space<Scalar> 
{
public:
L2_NEW_Space(bool serendipity= false);

L2_NEW_Space(MeshSharedPtr mesh, int p_init = 1,bool serendipity = false,Shapeset* shapeset = NULL);

virtual ~L2_NEW_Space();

virtual void set_shapeset(Shapeset* shapeset);

      /// Obtains an assembly list for the given element.
      virtual void get_element_assembly_list(Element* e, AsmList<Scalar>* al) const;


virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc);

/// Copy from Space instance 'space'
virtual void copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh);

protected:

virtual int get_edge_order(Element* e, int edge) const {
return H2D_MAKE_EDGE_ORDER(e->get_mode(), edge, this->edata[e->id].order);
}

bool serendipity; 

/// Common code for the constructors.
void init(Shapeset* shapeset, int p_init);


virtual void assign_bubble_dofs();

      /// Obtains an edge assembly list 
      virtual void get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al) const;


virtual void get_boundary_assembly_list_internal(Element* e, int ie, AsmList<Scalar>* al) const;
virtual void get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const;






friend class Space<Scalar>;

};
}
}









#endif
