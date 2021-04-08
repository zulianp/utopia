#ifndef UTOPIA_INTREPID2_SHELL_TOOLS_HPP
#define UTOPIA_INTREPID2_SHELL_TOOLS_HPP

#include <Intrepid2_CellTools.hpp>
#include <Kokkos_DynRankView.hpp>

#include "utopia_Views.hpp"

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar>
        class ShellTools {
        public:
            using ExecutionSpace = ::Kokkos::DefaultExecutionSpace;
            using DynRankView = ::Kokkos::DynRankView<Scalar>;
            using SizeType = ::Intrepid2::ordinal_type;

            // Jacobian : R^d-1 -> R^d
            // Jacobian^-1 : R^d -> R^d-1

            // Only affine elements
            template <class View>
            class ComputeJacobian {
            public:
                View cell_nodes;
                View jacobian;
                SizeType num_cells;
                SizeType num_qp;
                SizeType spatial_dim;

                UTOPIA_INLINE_FUNCTION ComputeJacobian(const View &cell_nodes, View &jacobian)
                    : cell_nodes(cell_nodes), jacobian(jacobian) {
                    num_cells = cell_nodes.extent(0);
                    spatial_dim = cell_nodes.extent(2);

                    assert(num_cells == jacobian.extent(0));
                    assert(spatial_dim == jacobian.extent(2));
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
                    for (SizeType d = 0; d < spatial_dim; ++d) {
                        for (SizeType n = 0; n < spatial_dim - 1; ++n) {
                            jacobian(cell, 0, d, n) = cell_nodes(cell, n + 1, d) - cell_nodes(cell, 0, d);
                        }
                    }

                    // FIXME, now it copies the first to the other q-points
                    for (SizeType qp = 1; qp < num_qp; ++qp) {
                        for (SizeType d = 0; d < spatial_dim; ++d) {
                            for (SizeType n = 0; n < spatial_dim - 1; ++n) {
                                jacobian(cell, qp, d, n) = jacobian(cell, 0, d, n);
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
                SizeType num_cells;
                SizeType num_qp;
                SizeType manifold_dim;
                SizeType spatial_dim;

                UTOPIA_INLINE_FUNCTION ComputeJacobianDeterminant(const View &jacobian, View &jacobian_det)
                    : jacobian(jacobian), jacobian_det(jacobian_det) {
                    num_cells = jacobian.extent(0);
                    num_qp = jacobian.extent(1);
                    manifold_dim = jacobian.extent(3);
                    spatial_dim = jacobian.extent(2);

                    assert(num_qp == jacobian_det.extent(1));
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
                    Scalar J = 0.0;

                    switch (manifold_dim) {
                        case 1: {
                            assert(spatial_dim == 2);
                            StaticMatrix<Scalar, 1, 2> mat;

                            for (SizeType d = 0; d < spatial_dim; ++d) {
                                mat(0, d) = jacobian(cell, 0, d, 0);
                            }

                            J = det(mat);
                            break;
                        }

                        case 2: {
                            assert(spatial_dim == 3);
                            StaticMatrix<Scalar, 3, 2> mat;
                            for (SizeType d1 = 0; d1 < spatial_dim; ++d1) {
                                for (SizeType d2 = 0; d2 < manifold_dim; ++d2) {
                                    mat(d1, d2) = jacobian(cell, 0, d1, d2);
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

                    for (SizeType qp = 0; qp < num_qp; ++qp) {
                        jacobian_det(cell, qp) = J;
                    }
                }
            };

            template <class View>
            class ComputeJacobianInverse {
            public:
                View jacobian;
                View jacobian_inv;
                SizeType num_cells;
                SizeType num_qp;
                SizeType manifold_dim;
                SizeType spatial_dim;

                UTOPIA_INLINE_FUNCTION ComputeJacobianInverse(const View &jacobian, View &jacobian_inv)
                    : jacobian(jacobian), jacobian_inv(jacobian_inv) {
                    num_cells = jacobian.extent(0);
                    num_qp = jacobian.extent(1);
                    spatial_dim = jacobian.extent(2);
                    manifold_dim = jacobian.extent(3);
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
                    switch (manifold_dim) {
                        case 1: {
                            assert(spatial_dim == 2);
                            StaticMatrix<Scalar, 2, 1> mat;
                            StaticMatrix<Scalar, 1, 2> mat_inv;

                            for (SizeType d = 0; d < spatial_dim; ++d) {
                                mat(d, 0) = jacobian(cell, 0, d, 0);
                            }

                            mat_inv = inv(mat);

                            for (SizeType d = 0; d < spatial_dim; ++d) {
                                jacobian_inv(cell, 0, 0, d) = mat_inv(0, d);
                            }

                            break;
                        }

                        case 2: {
                            assert(spatial_dim == 3);
                            StaticMatrix<Scalar, 3, 2> mat;
                            StaticMatrix<Scalar, 2, 3> mat_inv;

                            for (SizeType d1 = 0; d1 < spatial_dim; ++d1) {
                                for (SizeType d2 = 0; d2 < manifold_dim; ++d2) {
                                    mat(d1, d2) = jacobian(cell, 0, d1, d2);
                                }
                            }

                            mat_inv = inv(mat);

                            for (SizeType d1 = 0; d1 < manifold_dim; ++d1) {
                                for (SizeType d2 = 0; d2 < spatial_dim; ++d2) {
                                    jacobian_inv(cell, 0, d1, d2) = mat_inv(d1, d2);
                                }
                            }

                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }

                    for (SizeType qp = 0; qp < num_qp; ++qp) {
                        for (SizeType d1 = 0; d1 < manifold_dim; ++d1) {
                            for (SizeType d2 = 0; d2 < spatial_dim; ++d2) {
                                jacobian_inv(cell, qp, d1, d2) = jacobian_inv(cell, 0, d1, d2);
                            }
                        }
                    }
                }
            };

            template <class View>
            class AXPBY {
            public:
                Scalar a;
                View x;
                Scalar b;
                View y;

                SizeType components;

                UTOPIA_INLINE_FUNCTION AXPBY(const Scalar &a, const View &x, const Scalar &b, View &y)
                    : a(a), x(x), b(b), y(y) {
                    components = x.extent(1);
                    assert(components == y.extent(1));
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i) const {
                    for (SizeType c = 0; c < components; ++c) {
                        y(i, c) = a * x(i, c) + b * y(i, c);
                    }
                }
            };

            template <class View>
            class Scal {
            public:
                Scalar a;
                View x;

                SizeType components;

                UTOPIA_INLINE_FUNCTION Scal(const Scalar &a, View &x) : a(a), x(x) { components = x.extent(1); }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i) const {
                    for (SizeType c = 0; c < components; ++c) {
                        x(i, c) *= a;
                    }
                }
            };

            template <class View>
            class TransformGradientToPhysicalSpace {
            public:
                View jacobian_inv;
                View grad_ref;
                View grad_physical;

                SizeType manifold_dim;

                TransformGradientToPhysicalSpace(const View &jacobian_inv, const View &grad_ref, View &grad_physical)
                    : jacobian_inv(jacobian_inv), grad_ref(grad_ref), grad_physical(grad_physical) {
                    manifold_dim = jacobian_inv.extent(2);
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell,
                                                       const SizeType node,
                                                       const SizeType qp,
                                                       const SizeType d_physical) const {
                    for (SizeType d_ref = 0; d_ref < manifold_dim; ++d_ref) {
                        grad_physical(cell, node, qp, d_physical) +=
                            jacobian_inv(cell, qp, d_ref, d_physical) * grad_ref(node, qp, d_ref);
                    }
                }
            };

            static void allocate_jacobian(const SizeType n_cells,
                                          const SizeType spatial_dim,
                                          const SizeType n_qp,
                                          DynRankView &jacobian) {
                jacobian = DynRankView("jacobian", n_cells, n_qp, spatial_dim, spatial_dim - 1);
            }

            static void allocate_jacobian_inverse(const SizeType n_cells,
                                                  const SizeType spatial_dim,
                                                  const SizeType n_qp,
                                                  DynRankView &jacobian_inv) {
                jacobian_inv = DynRankView("jacobian_inv", n_cells, n_qp, spatial_dim - 1, spatial_dim);
            }

            static void allocate_measure(const SizeType n_cells, const SizeType n_qp, DynRankView &measure) {
                measure = DynRankView("measure", n_cells, n_qp);
            }

            static void transform_gradient_to_physical_space(const DynRankView &jacobian_inv,
                                                             const DynRankView &grad_ref,
                                                             DynRankView &grad_physical) {
                const SizeType num_cells = jacobian_inv.extent(0);
                const SizeType n_qp = jacobian_inv.extent(1);
                const SizeType num_nodes = grad_ref.extent(0);
                const SizeType spatial_dim = jacobian_inv.extent(1);

                TransformGradientToPhysicalSpace<DynRankView> transform(jacobian_inv, grad_ref, grad_physical);

                Kokkos::parallel_for("ShellTools::transform_gradient_to_physical_space",
                                     Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>(
                                         {0, 0, 0, 0}, {num_cells, num_nodes, n_qp, spatial_dim}),
                                     transform);
            }

            static void cell_geometry(const DynRankView &cell_nodes,
                                      DynRankView &jacobian,
                                      DynRankView &jacobian_inv,
                                      DynRankView &measure) {
                ComputeJacobian<DynRankView> compute_J(cell_nodes, jacobian);
                ComputeJacobianInverse<DynRankView> compute_J_inv(jacobian, jacobian_inv);
                ComputeJacobianDeterminant<DynRankView> compute_J_det(jacobian, measure);

                const SizeType spatial_dim = cell_nodes.extent(2);
                Scal<DynRankView> scal((spatial_dim == 2) ? 0.5 : 1.0 / 6.0, measure);

                SizeType num_cells = jacobian.extent(0);
                Kokkos::parallel_for(
                    "ShellTools::cell_geometry", num_cells, KOKKOS_LAMBDA(const SizeType &cell) {
                        compute_J(cell);
                        compute_J_det(cell);  // This depends on compute_J(cell)
                        compute_J_inv(cell);
                        scal(cell);
                    });
            }
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_SHELL_TOOLS_HPP
