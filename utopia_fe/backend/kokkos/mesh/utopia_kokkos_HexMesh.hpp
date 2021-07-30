#ifndef UTOPIA_KOKKOS_HEX_MESH_HPP
#define UTOPIA_KOKKOS_HEX_MESH_HPP

#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace kokkos {

        template <typename Scalar, class ExecutionSpace>
        class HexMesh {
        public:
            void init(const int nx, const int ny, const int nz) {
                nodeCoord = Kokkos::DynRankView<Scalar, ExecutionSpace>("nodeCoord", numNodes, spaceDim);
                nodeOnBoundary = Kokkos::DynRankView<int, ExecutionSpace>("nodeOnBoundary", numNodes);

                int inode = 0;
                for (int k = 0; k < nz + 1; k++) {
                    for (int j = 0; j < ny + 1; j++) {
                        for (int i = 0; i < nx + 1; i++) {
                            nodeCoord(inode, 0) = leftX + (double)i * hx;
                            nodeCoord(inode, 1) = leftY + (double)j * hy;
                            nodeCoord(inode, 2) = leftZ + (double)k * hz;
                            if (k == 0 || j == 0 || i == 0 || k == nz || j == ny || i == nx) {
                                // FIXME USE correct boundary tag
                                nodeOnBoundary(inode) = 1;
                            } else {
                                nodeOnBoundary(inode) = 0;
                            }
                            inode++;
                        }
                    }
                }

                elemToNode = Kokkos::DynRankView<int, DeviceSpaceType>("elemToNode", numElems, numNodesPerElem);

                int ielem = 0;
                for (int k = 0; k < nz; k++) {
                    for (int j = 0; j < ny; j++) {
                        for (int i = 0; i < nx; i++) {
                            elemToNode(ielem, 0) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i;
                            elemToNode(ielem, 1) = (ny + 1) * (nx + 1) * k + (nx + 1) * j + i + 1;
                            elemToNode(ielem, 2) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i + 1;
                            elemToNode(ielem, 3) = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + i;
                            elemToNode(ielem, 4) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i;
                            elemToNode(ielem, 5) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + i + 1;
                            elemToNode(ielem, 6) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i + 1;
                            elemToNode(ielem, 7) = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + i;
                            ielem++;
                        }
                    }
                }
            }

            Kokkos::DynRankView<Scalar, ExecutionSpace> nodeCoord;
            Kokkos::DynRankView<int, ExecutionSpace> elemToNode;
            Kokkos::DynRankView<int, ExecutionSpace> nodeOnBoundary;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_HEX_MESH_HPP
