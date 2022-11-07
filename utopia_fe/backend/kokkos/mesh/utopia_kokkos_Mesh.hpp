#ifndef UTOPIA_KOKKOS_MESH_HPP
#define UTOPIA_KOKKOS_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Input.hpp"
#include "utopia_kokkos_Cube.hpp"

#include <Kokkos_DynRankView.hpp>

#include <vector>

namespace utopia {
    namespace kokkos {

        template <typename Scalar, class ExecutionSpace_ = ::Kokkos::DefaultExecutionSpace, typename IntType = int>
        class Mesh : public Configurable {
        public:
#ifdef UTOPIA_WITH_MPI
            using Comm = utopia::MPICommunicator;
#else
            using Comm = utopia::SelfCommunicator;
#endif
            using ExecutionSpace = ExecutionSpace_;

            void read(Input &in) override {
                std::string type = "cube";
                in.get("type", type);

                if (type == "cube") {
                    Cube<Mesh> cube;
                    cube.read(in);
                    cube.build(*this);
                } else {
                    Utopia::Abort();
                }
            }

            //
            class Part {
            public:
                using ExecutionSpace = ExecutionSpace_;

                void set_elem_dims(const int n, const int nnodexelement) {
                    elem_to_node =
                        Kokkos::DynRankView<IntType, ExecutionSpace>("Mesh::Part::elem_to_node", n, nnodexelement);
                }

                // Local indexing
                Kokkos::DynRankView<IntType, ExecutionSpace> elem_to_node;
                int part_id{-1};
                int elem_type{-1};
                std::string name;
            };

            void set_point_dims(int n, int dims) {
                points = Kokkos::DynRankView<Scalar, ExecutionSpace>("Mesh::points", n, dims);
            }

            void set_n_parts(int n) { parts.resize(n); }

            // Mesh(const Comm &comm)
            // : comm(comm) {}

            Comm comm;
            Kokkos::DynRankView<Scalar, ExecutionSpace> points;
            std::vector<Part> parts;
        };
    }  // namespace kokkos

}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MESH_HPP
