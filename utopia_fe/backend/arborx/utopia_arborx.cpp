#include <ArborX.hpp>
#include <ArborX_DistributedTree.hpp>
#include <ArborX_Ray.hpp>
#include <ArborX_Version.hpp>

// https://github.com/arborx/ArborX/blob/master/benchmarks/distributed_tree_driver/distributed_tree_driver.cpp

namespace utopia {
    namespace arborx {

        void test() {
            using HostExecutionSpace = ::Kokkos::DefaultHostExecutionSpace;
            using ExecutionSpace = Kokkos::Serial;
            // FIXME
            using DeviceType = Kokkos::Serial;
            using MemorySpace = DeviceType::memory_space;

            auto comm = MPI_COMM_WORLD;

            ::Kokkos::View<::ArborX::Box *, DeviceType> bounding_boxes("Testing::bounding_boxes");
            ::ArborX::DistributedTree<MemorySpace> distributed_tree(comm, ExecutionSpace{}, bounding_boxes);
        }
    }  // namespace arborx
}  // namespace utopia
