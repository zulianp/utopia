#ifndef UTOPIA_PETSC_MAKE_ELEM_HPP
#define UTOPIA_PETSC_MAKE_ELEM_HPP

#include <cassert>
#include "utopia_ElemTraits.hpp"
#include "utopia_Tri3.hpp"

namespace utopia {

    template <class Space, class Elem, bool IsSimplex = is_simplex<Elem>::value>
    class MakeElem
    // {};

    // template<class Space, class Elem>
    // class MakeElem<Space, Elem, false>
    {
    public:
        using SizeType = typename Space::SizeType;

        inline static void apply(const Space &space, const SizeType &idx, Elem &e) {
            const auto &mesh = space.mesh();
            // mesh.elem(idx, e.univar_elem());

            e.idx(idx);
            typename Space::Mesh::Point translation, cell_size;
            mesh.cell_point(idx, translation);
            mesh.cell_size(idx, cell_size);
            e.set(translation, cell_size);
        }
    };

    template <class Space, typename Scalar, int PhysicalDim, int NVar>
    class MakeElem<Space, MultiVariateElem<Tri3<Scalar, PhysicalDim>, NVar>, true> {
    public:
        using SizeType = typename Space::SizeType;
        using Point = typename Space::Point;
        using Elem = utopia::MultiVariateElem<Tri3<Scalar, PhysicalDim>, NVar>;

        inline static void apply(const Space &space, const SizeType &idx, Elem &e) {
            const auto &mesh = space.mesh();
            // mesh.elem(idx, e.univar_elem());
            e.univar_elem().idx(idx);

            typename Space::Mesh::NodeIndex nodes;
            mesh.nodes(idx, nodes);

            const SizeType n = nodes.size();

            Point p0, p1, p2;

            if (n == 3) {
                mesh.point(nodes[0], p0);
                mesh.point(nodes[1], p1);
                mesh.point(nodes[2], p2);

                e.set(p0, p1, p2);
            } else {
                assert(false);
                // Point p4;

                // e.init()
            }
        }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_MAKE_ELEM_HPP
