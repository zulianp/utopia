#ifndef UTOPIA_KOKKOS_HEX_MESH_HPP
#define UTOPIA_KOKKOS_HEX_MESH_HPP

#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace kokkos {

        template <typename Scalar, class ExecutionSpace>
        class HexMesh {
        public:
            void init(const int nx, const int ny, const int nz) {
                node_coord = Kokkos::DynRankView<Scalar, ExecutionSpace>("node_coord", numNodes, spaceDim);
                node_on_boundary = Kokkos::DynRankView<int, ExecutionSpace>("node_on_boundary", numNodes);

                int inode = 0;
                for (int k = 0; k < nz + 1; k++) {
                    for (int j = 0; j < ny + 1; j++) {
                        for (int i = 0; i < nx + 1; i++) {
                            node_coord(inode, 0) = leftX + (double)i * hx;
                            node_coord(inode, 1) = leftY + (double)j * hy;
                            node_coord(inode, 2) = leftZ + (double)k * hz;
                            if (k == 0 || j == 0 || i == 0 || k == nz || j == ny || i == nx) {
                                // FIXME USE correct boundary tag
                                node_on_boundary(inode) = 1;
                            } else {
                                node_on_boundary(inode) = 0;
                            }
                            inode++;
                        }
                    }
                }

                elem_to_node = Kokkos::DynRankView<int, DeviceSpaceType>("elem_to_node", numElems, numNodesPerElem);

                int ielem = 0;
                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            elem_to_node(ielem, 0) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i;
                            elem_to_node(ielem, 1) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i + 1;
                            elem_to_node(ielem, 2) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i + 1;
                            elem_to_node(ielem, 3) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i;
                            elem_to_node(ielem, 4) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i;
                            elem_to_node(ielem, 5) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i + 1;
                            elem_to_node(ielem, 6) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i + 1;
                            elem_to_node(ielem, 7) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i;
                            ielem++;
                        }
                    }
                }
            }

            Kokkos::DynRankView<Scalar, ExecutionSpace> node_coord;
            Kokkos::DynRankView<int, ExecutionSpace> elem_to_node;
            Kokkos::DynRankView<int, ExecutionSpace> node_on_boundary;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_HEX_MESH_HPP
