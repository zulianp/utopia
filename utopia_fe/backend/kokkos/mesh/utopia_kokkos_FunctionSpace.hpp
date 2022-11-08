#ifndef UTOPIA_KOKKOS_FUNCTIONSPACE_HPP
#define UTOPIA_KOKKOS_FUNCTIONSPACE_HPP

#include "utopia_kokkos_Mesh.hpp"

namespace utopia {
    namespace kokkos {
        template <class ExecutionSpace_, typename IntType>
        class DofMap {
        public:
            using ExecutionSpace = ExecutionSpace_;
            // Map mesh local-indexing to dof global indexing
            Kokkos::DynRankView<IntType, ExecutionSpace> dof_map;

            // In case no dof_map we pass this!
            class IdentityMap {
            public:
                UTOPIA_INLINE_FUNCTION IntType operator()(const IntType idx) const { return idx; }
            };
        };

        template <typename Scalar, class ExecutionSpace_ = ::Kokkos::DefaultExecutionSpace, typename IntType = int>
        class FunctionSpace : public Configurable {
        public:
            using Mesh = utopia::kokkos::Mesh<Scalar, ExecutionSpace_, IntType>;
            using DofMap = utopia::kokkos::DofMap<ExecutionSpace_, IntType>;
            using Comm = typename Mesh::Comm;
            using ExecutionSpace = ExecutionSpace_;
            using Part = typename Mesh::Part;

            void read(Input &in) override {
                in.get("mesh", mesh);
                parts = mesh.parts;
            }

            FunctionSpace(const Comm &comm) : mesh(comm) {}

            Mesh mesh;
            std::vector<Part> parts;

            void create_matrix() {}
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FUNCTIONSPACE_HPP
