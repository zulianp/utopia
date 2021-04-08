#ifndef UTOPIA_INTREPID2_SHELL_TOOLS_HPP
#define UTOPIA_INTREPID2_SHELL_TOOLS_HPP

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar>
        class ShellTools {
        public:
            // using ExecutionSpace = ::Kokkos::Serial;
            using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;
            using DynRankView = ::Kokkos::DynRankView<Scalar>;
            using Cubature = ::Intrepid2::Cubature<ExecutionSpace, Scalar, Scalar>;
            using CubaturePtr = ::Teuchos::RCP<Cubature>;
            using CellTopology = ::shards::CellTopology;
            // using Tet = ::Intrepid2::Basis_HGRAD_TET_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Tri = ::Intrepid2::Basis_HGRAD_TRI_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using Line = ::Intrepid2::Basis_HGRAD_LINE_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            // using Hex = ::Intrepid2::Basis_HGRAD_HEX_C1_FEM<ExecutionSpace, Scalar, Scalar>;
            using CellTools = ::Intrepid2::CellTools<ExecutionSpace>;
            using FunctionSpaceTools = ::Intrepid2::FunctionSpaceTools<ExecutionSpace>;
            using SizeType = ::Intrepid2::ordinal_type;

            // Functors

            // Only affine elements

            template <class View>
            class ComputeJacobian {
            public:
                View cell_nodes;
                View jacobian;
                int num_cells;
                int num_nodes;
                int spatial_dim;
                int num_qp;

                UTOPIA_INLINE_FUNCTION ComputeJacobian(const View &cell_nodes, View &jacobian)
                    : cell_nodes(cell_nodes), jacobian(jacobian) {
                    num_cells = cell_nodes.extent(0);
                    num_nodes = cell_nodes.extent(1);
                    spatial_dim = cell_nodes.extent(2);
                    num_qp = jacobian.exxtent(3);

                    assert(num_nodes == spatial_dim);

                    assert(num_cells == jacobian.extent(0));
                    assert(num_nodes - 1 == jacobian.extent(1));
                    assert(spatial_dim == jacobian.extent(2));
                }

                UTOPIA_INLINE_FUNCTION void operator()(const int cell) const {
                    for (int n = 0; n < num_nodes - 1; ++n) {
                        for (int d = 0; d < spatial_dim; ++d) {
                            jacobian(cell, n, d, 0) = cell_nodes(cell, n + 1, d) - cell_nodes(cell, 0, d);
                        }
                    }

                    // copy to other q-points (fix for warped elements)
                    for (int n = 0; n < num_nodes - 1; ++n) {
                        for (int d = 0; d < spatial_dim; ++d) {
                            for (int qp = 1; qp < num_qp; ++qp) {
                                jacobian(cell, n, d, qp) = jacobian(cell, n, d, 0);
                            }
                        }
                    }
                }
            };

            template <class View>
            class ComputeJacobianDeterminant {
            public:
                View jacobian;
                View jacobian_det;
                int num_cells;
                int manifold_dim;
                int spatial_dim;
                int num_qp;

                UTOPIA_INLINE_FUNCTION ComputeJacobianDeterminant(const View &jacobian, View &jacobian_det)
                    : jacobian(jacobian), jacobian_det(jacobian_det) {
                    num_cells = jacobian.extent(0);
                    manifold_dim = jacobian.extent(1);
                    spatial_dim = jacobian.extent(2);
                    num_qp = jacobian.extent(3);

                    assert(num_qp == jacobian_det.extent(1));
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const int cell) const {
                    Scalar J = 0.0;

                    switch (manifold_dim) {
                        case 1: {
                            assert(spatial_dim == 3);
                            StaticMatrix<Scalar, 1, 2> mat;

                            for (int d = 0; d < spatial_dim; ++d) {
                                mat(0, d) = jacobian(cell, 0, d, 0);
                            }

                            J = det(mat);
                            break;
                        }

                        case 2: {
                            assert(spatial_dim == 3);
                            StaticMatrix<Scalar, 3, 2> mat;
                            for (int d1 = 0; d1 < manifold_dim; ++d1) {
                                for (int d2 = 0; d2 < spatial_dim; ++d2) {
                                    mat(d1, d2) = jacobian(cell, d1, d2, 0);
                                }
                            }

                            J = det(mat);
                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }

                    for (int qp = 0; qp < num_qp; ++qp) {
                        jacobian_det(cell, qp) = J;
                    }
                }
            };

            void cell_geometry(const DynRankView &cell_nodes, DynRankView &jacobian, DynRankView &jacobian_det) {
                ComputeJacobianDeterminant<DynRankView> compute_J(cell_nodes, jacobian);
                ComputeJacobianDeterminant<DynRankView> compute_J_det(jacobian, jacobian_det);

                int num_cells = jacobian.extent(0);
                Kokkos::parallel_for(
                    "ShellTools::cell_geometry", num_cells, KOKKOS_LAMBDA(const int &cell) {
                        compute_J(cell);
                        compute_J_det(cell);  // This depends on compute_J(cell)
                    });
            }
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_SHELL_TOOLS_HPP
