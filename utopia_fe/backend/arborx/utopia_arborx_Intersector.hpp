#ifndef UTOPIA_ARBORX_INTERSECTOR_HPP
#define UTOPIA_ARBORX_INTERSECTOR_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_fe_base.hpp"

#include <ArborX_DistributedTree.hpp>
#include <Kokkos_Core.hpp>

namespace utopia {
    namespace arborx {

        template <typename DeviceType>
        struct BoxSearches {
            ::Kokkos::View<::ArborX::Box *, DeviceType> boxes;
        };

    }  // namespace arborx
}  // namespace utopia

template <typename DeviceType>
struct ArborX::AccessTraits<utopia::arborx::BoxSearches<DeviceType>, ArborX::PredicatesTag> {
    using memory_space = typename DeviceType::memory_space;
    static KOKKOS_FUNCTION std::size_t size(utopia::arborx::BoxSearches<DeviceType> const &pred) {
        return pred.boxes.extent(0);
    }
    static KOKKOS_FUNCTION auto get(utopia::arborx::BoxSearches<DeviceType> const &pred, std::size_t i) {
        return ArborX::intersects(pred.boxes(i));
    }
};

namespace utopia {
    namespace arborx {

        template <class ElemView>
        class DetectIntersectionsFromElementArray {
        public:
            using HostExecutionSpace = ::Kokkos::DefaultHostExecutionSpace;
            using ExecutionSpace = Kokkos::Serial;
            // FIXME
            using DeviceType = Kokkos::Serial;
            using MemorySpace = DeviceType::memory_space;
            using ArborXScalar = float;
            using PairIndexRank = Kokkos::pair<int, int>;
            using BoxSearches = utopia::arborx::BoxSearches<DeviceType>;

            DetectIntersectionsFromElementArray(MPI_Comm comm)
                : comm_(comm), offsets("Testing::offsets", 0), values("Testing::values", 0) {}

            void create_bounding_boxes(const ElemView &cell_nodes,
                                       ::Kokkos::View<::ArborX::Box *, DeviceType> &bounding_boxes) {
                auto n_elems = cell_nodes.extent(0);
                int n_points = cell_nodes.extent(1);
                int dim = cell_nodes.extent(2);

                Kokkos::parallel_for(
                    "DetectIntersectionsFromElementArray::construct_bounding_boxes",
                    Kokkos::RangePolicy<ExecutionSpace>(0, n_elems),
                    KOKKOS_LAMBDA(int i) {
                        ArborX::Box box;
                        for (int k = 0; k < n_points; ++k) {
                            for (int d = 0; d < dim; ++d) {
                                box.minCorner()[d] = device::min(box.minCorner()[d], ArborXScalar(cell_nodes(i, k, d)));
                                box.maxCorner()[d] = device::max(box.maxCorner()[d], ArborXScalar(cell_nodes(i, k, d)));
                            }
                        }

                        bounding_boxes(i) = box;
                    });
            }

            bool detect(const ElemView &from, const ElemView &to) {
                ::Kokkos::View<::ArborX::Box *, DeviceType> bounding_boxes_from(
                    ::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                         "DetectIntersectionsFromElementArray::bounding_boxes"),
                    from.extent(0));

                create_bounding_boxes(from, bounding_boxes_from);

                ::Kokkos::View<::ArborX::Box *, DeviceType> bounding_boxes_to(
                    ::Kokkos::view_alloc(::Kokkos::WithoutInitializing,
                                         "DetectIntersectionsFromElementArray::bounding_boxes"),
                    to.extent(0));

                create_bounding_boxes(to, bounding_boxes_to);

                ::ArborX::DistributedTree<MemorySpace> distributed_tree(comm_, ExecutionSpace{}, bounding_boxes_from);

                Kokkos::View<int *, DeviceType> offsets("Testing::offsets", 0);
                Kokkos::View<PairIndexRank *, DeviceType> values("Testing::values", 0);

                BoxSearches box_searches;
                box_searches.boxes = bounding_boxes_to;
                distributed_tree.query(ExecutionSpace{}, box_searches, values, offsets);

                return false;
            }

            void describe(std::ostream &os) const {
                os << "values.extent(0) = " << values.extent(0) << "\n";

                for (int i = 0; i < values.extent(0); ++i) {
                    os << values(i).first << " -> " << values(i).second << "\n";
                }

                os << "offsets.extent(0) = " << offsets.extent(0) << "\n";

                for (int i = 0; i < offsets.extent(0); ++i) {
                    os << offsets(i) << "\n";
                }
            }

        private:
            MPI_Comm comm_;

            Kokkos::View<int *, DeviceType> offsets;
            Kokkos::View<PairIndexRank *, DeviceType> values;
        };
    }  // namespace arborx
}  // namespace utopia

#endif  // UTOPIA_ARBORX_INTERSECTOR_HPP
