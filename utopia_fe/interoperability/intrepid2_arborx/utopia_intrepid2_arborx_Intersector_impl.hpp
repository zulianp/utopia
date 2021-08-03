#ifndef UTOPIA_INTREPID2_ARBORX_INTERSECTOR_IMPL_HPP
#define UTOPIA_INTREPID2_ARBORX_INTERSECTOR_IMPL_HPP

#include "utopia_arborx_Intersector.hpp"
#include "utopia_intrepid2_arborx_Intersector.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {

    template <class CellPoints>
    class Intrepid2ArborXIntersector<CellPoints>::Impl {
    public:
    };

    template <class CellPoints>
    Intrepid2ArborXIntersector<CellPoints>::Intrepid2ArborXIntersector() : impl_(utopia::make_unique<Impl>()) {}

    template <class CellPoints>
    Intrepid2ArborXIntersector<CellPoints>::~Intrepid2ArborXIntersector() = default;

    template <class CellPoints>
    bool Intrepid2ArborXIntersector<CellPoints>::detect(const CellPoints &points1, const CellPoints &points2) {
        arborx::DetectIntersectionsFromElementArray<CellPoints> isect(MPI_COMM_WORLD);
        if (!isect.detect(points1, points2)) {
            return false;
        }

        return true;
    }

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_ARBORX_INTERSECTOR_IMPL_HPP
