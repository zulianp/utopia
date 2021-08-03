#ifndef UTOPIA_INTREPID2_SHELL_TOOLS_HPP
#define UTOPIA_INTREPID2_SHELL_TOOLS_HPP

#include "utopia_intrepid2_Base.hpp"

#include "utopia_Views.hpp"

#include <Intrepid2_CellTools.hpp>
#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar>
        class ShellTools {
        public:
            using HostExecutionSpace = utopia::intrepid2::HostExecutionSpace;
            using ExecutionSpace = utopia::intrepid2::ExecutionSpace;
            using DynRankView = utopia::intrepid2::ViewDevice<Scalar>;
            using IntView = utopia::intrepid2::IntViewDevice;

            using SizeType = ::Intrepid2::ordinal_type;
            using CellRange = Kokkos::RangePolicy<ExecutionSpace>;
            using CellQPRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;

            // Jacobian : R^d-1 -> R^d
            // Jacobian^-1 : R^d -> R^d-1

            // Only simplexes
            // template <class View>
            // class SimplexComputeJacobian {
            // public:
            //     const View cell_nodes;
            //     View jacobian;
            //     const SizeType n_cells;
            //     const SizeType n_quad_points;
            //     const SizeType spatial_dim;

            //     UTOPIA_INLINE_FUNCTION SimplexComputeJacobian(const View &cell_nodes, View &jacobian)
            //         : cell_nodes(cell_nodes),
            //           jacobian(jacobian),
            //           n_cells(cell_nodes.extent(0)),
            //           n_quad_points(jacobian.extent(1)),
            //           spatial_dim(cell_nodes.extent(2)) {
            //         assert(n_cells == jacobian.extent(0));
            //         assert(spatial_dim == jacobian.extent(2));
            //     }

            //     UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
            //         for (SizeType d = 0; d < spatial_dim; ++d) {
            //             for (SizeType n = 0; n < spatial_dim - 1; ++n) {
            //                 jacobian(cell, 0, d, n) = cell_nodes(cell, n + 1, d) - cell_nodes(cell, 0, d);
            //             }
            //         }

            //         // FIXME, now it copies the first to the other q-points
            //         for (SizeType qp = 1; qp < n_quad_points; ++qp) {
            //             for (SizeType d = 0; d < spatial_dim; ++d) {
            //                 for (SizeType n = 0; n < spatial_dim - 1; ++n) {
            //                     jacobian(cell, qp, d, n) = jacobian(cell, 0, d, n);
            //                 }
            //             }
            //         }
            //     }
            // };

            template <class View>
            class IsoParametricComputeJacobian {
            public:
                const View cell_nodes;
                const View ref_grad;
                View jacobian;
                const SizeType n_cells;
                const SizeType spatial_dim;
                const SizeType num_nodes;

                UTOPIA_INLINE_FUNCTION IsoParametricComputeJacobian(const View &cell_nodes,
                                                                    const View &ref_grad,
                                                                    View &jacobian)
                    : cell_nodes(cell_nodes),
                      ref_grad(ref_grad),
                      jacobian(jacobian),
                      n_cells(cell_nodes.extent(0)),
                      spatial_dim(cell_nodes.extent(2)),
                      num_nodes(cell_nodes.extent(1)) {
                    assert(n_cells == SizeType(jacobian.extent(0)));
                    assert(spatial_dim == SizeType(jacobian.extent(2)));
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell, const SizeType qp) const {
                    for (SizeType sdim = 0; sdim < spatial_dim; ++sdim) {
                        for (SizeType mdim = 0; mdim < spatial_dim - 1; ++mdim) {
                            jacobian(cell, qp, sdim, mdim) = 0.0;
                        }
                    }

                    for (SizeType i = 0; i < num_nodes; ++i) {
                        for (SizeType sdim = 0; sdim < spatial_dim; ++sdim) {
                            for (SizeType mdim = 0; mdim < spatial_dim - 1; ++mdim) {
                                jacobian(cell, qp, sdim, mdim) += cell_nodes(cell, i, sdim) * ref_grad(i, qp, mdim);
                            }
                        }
                    }
                }
            };

            template <class View>
            class ComputeJacobianDeterminant {
            public:
                const View jacobian;
                View jacobian_det;
                const SizeType n_cells;
                const SizeType n_quad_points;
                const SizeType manifold_dim;
                const SizeType spatial_dim;

                UTOPIA_INLINE_FUNCTION ComputeJacobianDeterminant(const View &jacobian, View &jacobian_det)
                    : jacobian(jacobian),
                      jacobian_det(jacobian_det),
                      n_cells(jacobian.extent(0)),
                      n_quad_points(jacobian.extent(1)),
                      manifold_dim(jacobian.extent(3)),
                      spatial_dim(jacobian.extent(2)) {
                    assert(n_quad_points == SizeType(jacobian_det.extent(1)));
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
                    (*this)(cell, 0);

                    const Scalar J = jacobian_det(cell, 0);
                    for (SizeType qp = 0; qp < n_quad_points; ++qp) {
                        jacobian_det(cell, qp) = J;
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell, const SizeType qp) const {
                    Scalar J = 0.0;

                    switch (manifold_dim) {
                        case 1: {
                            J = 0.0;
                            for (SizeType d = 0; d < spatial_dim; ++d) {
                                auto x = jacobian(cell, qp, d, 0);
                                J += x * x;
                            }

                            J = device::sqrt(J);
                            assert(J != 0);
                            break;
                        }

                        case 2: {
                            assert(spatial_dim == 3);
                            StaticMatrix<Scalar, 3, 2> mat;

                            for (SizeType d1 = 0; d1 < 3; ++d1) {
                                for (SizeType d2 = 0; d2 < 2; ++d2) {
                                    mat(d1, d2) = jacobian(cell, qp, d1, d2);
                                }
                            }

                            J = det(mat);
                            assert(J != 0);

                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }

                    jacobian_det(cell, qp) = J;
                }
            };

            template <class View>
            class ComputeJacobianInverse {
            public:
                const View jacobian;
                View jacobian_inv;
                const SizeType n_cells;
                const SizeType n_quad_points;
                const SizeType spatial_dim;
                const SizeType manifold_dim;

                UTOPIA_INLINE_FUNCTION ComputeJacobianInverse(const View &jacobian, View &jacobian_inv)
                    : jacobian(jacobian),
                      jacobian_inv(jacobian_inv),
                      n_cells(jacobian.extent(0)),
                      n_quad_points(jacobian.extent(1)),
                      spatial_dim(jacobian.extent(2)),
                      manifold_dim(jacobian.extent(3)) {}

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell) const {
                    (*this)(cell, 0);

                    for (SizeType qp = 0; qp < n_quad_points; ++qp) {
                        for (SizeType d1 = 0; d1 < manifold_dim; ++d1) {
                            for (SizeType d2 = 0; d2 < spatial_dim; ++d2) {
                                jacobian_inv(cell, qp, d1, d2) = jacobian_inv(cell, 0, d1, d2);
                            }
                        }
                    }
                }

                template <int SpatialDim, int ManifoldDim>
                UTOPIA_INLINE_FUNCTION void compute_inv(const SizeType cell, const SizeType qp) const {
                    StaticMatrix<Scalar, SpatialDim, ManifoldDim> mat;
                    StaticMatrix<Scalar, ManifoldDim, SpatialDim> mat_inv;

                    for (SizeType d1 = 0; d1 < SpatialDim; ++d1) {
                        for (SizeType d2 = 0; d2 < ManifoldDim; ++d2) {
                            mat(d1, d2) = jacobian(cell, qp, d1, d2);
                        }
                    }

                    mat_inv = inv(mat);

                    assert(approxeq(det(mat_inv), 1. / det(mat), 1e-8));

                    for (SizeType d1 = 0; d1 < ManifoldDim; ++d1) {
                        for (SizeType d2 = 0; d2 < SpatialDim; ++d2) {
                            jacobian_inv(cell, qp, d1, d2) = mat_inv(d1, d2);
                        }
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell, const SizeType qp) const {
                    switch (manifold_dim) {
                        case 1: {
                            if (spatial_dim == 2) {
                                this->compute_inv<2, 1>(cell, qp);
                            } else if (spatial_dim == 3) {
                                this->compute_inv<3, 1>(cell, qp);
                            }

                            break;
                        }

                        case 2: {
                            assert(spatial_dim == 3);
                            this->compute_inv<3, 2>(cell, qp);
                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }
                }
            };

            template <class View>
            class AXPBY {
            public:
                const Scalar a;
                const View x;
                const Scalar b;
                View y;

                const SizeType components;

                UTOPIA_INLINE_FUNCTION AXPBY(const Scalar &a, const View &x, const Scalar &b, View &y)
                    : a(a), x(x), b(b), y(y), components(x.extent(1)) {
                    assert(components == y.extent(1));
                }

                // Only affine elements
                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i) const {
                    for (SizeType c = 0; c < components; ++c) {
                        y(i, c) = a * x(i, c) + b * y(i, c);
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i, const SizeType qp) const {
                    y(i, qp) = a * x(i, qp) + b * y(i, qp);
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

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i, const SizeType qp) const { x(i, qp) *= a; }
            };

            template <class View>
            class ScalJac {
            public:
                const Scalar a;
                View J;

                const SizeType rows, cols;

                UTOPIA_INLINE_FUNCTION ScalJac(const Scalar &a, View &J)
                    : a(a), J(J), rows(J.extent(2)), cols(J.extent(3)) {}

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType i, const SizeType qp) const {
                    for (SizeType r = 0; r < rows; ++r) {
                        for (SizeType c = 0; c < cols; ++c) {
                            J(i, qp, r, c) *= a;
                        }
                    }
                }
            };

            template <class View>
            class CheckGrad {
            public:
                CheckGrad(const View &grad_physical)
                    : grad_physical(grad_physical), spatial_dim(grad_physical.extent(3)) {}

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell,
                                                       const SizeType node,
                                                       const SizeType qp) const

                {
                    bool has_nnz = false;
                    for (SizeType d = 0; d < spatial_dim; ++d) {
                        if (device::abs(grad_physical(cell, node, qp, d)) > 0) {
                            has_nnz = true;
                        }
                    }

                    assert(has_nnz);
                }

                View grad_physical;
                SizeType spatial_dim;
            };

            static void check_grad(const DynRankView &grad_physical) {
                CheckGrad<DynRankView> check(grad_physical);

                const SizeType n_cells = grad_physical.extent(0);
                const SizeType num_nodes = grad_physical.extent(1);
                const SizeType n_qp = grad_physical.extent(2);

                Kokkos::parallel_for(
                    "ShellTools::check_grad",
                    Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {n_cells, num_nodes, n_qp}),
                    check);
            }

            template <class View>
            class TransformGradientToPhysicalSpace {
            public:
                View jacobian_inv;
                View ref_grad;
                View grad_physical;

                SizeType manifold_dim;

                TransformGradientToPhysicalSpace(const View &jacobian_inv, const View &ref_grad, View &grad_physical)
                    : jacobian_inv(jacobian_inv), ref_grad(ref_grad), grad_physical(grad_physical) {
                    manifold_dim = jacobian_inv.extent(2);
                }

                UTOPIA_INLINE_FUNCTION void operator()(const SizeType cell,
                                                       const SizeType node,
                                                       const SizeType qp,
                                                       const SizeType d_physical) const {
                    for (SizeType d_ref = 0; d_ref < manifold_dim; ++d_ref) {
                        grad_physical(cell, node, qp, d_physical) +=
                            jacobian_inv(cell, qp, d_ref, d_physical) * ref_grad(node, qp, d_ref);
                    }
                }
            };

            static void allocate_jacobian(const SizeType n_cells,
                                          const SizeType manifold_dim,
                                          const SizeType spatial_dim,
                                          const SizeType n_qp,
                                          DynRankView &jacobian) {
                jacobian = DynRankView("jacobian", n_cells, n_qp, spatial_dim, manifold_dim);
            }

            static void allocate_jacobian_inverse(const SizeType n_cells,
                                                  const SizeType manifold_dim,
                                                  const SizeType spatial_dim,
                                                  const SizeType n_qp,
                                                  DynRankView &jacobian_inv) {
                jacobian_inv = DynRankView("jacobian_inv", n_cells, n_qp, manifold_dim, spatial_dim);
            }

            static void allocate_measure(const SizeType n_cells, const SizeType n_qp, DynRankView &measure) {
                measure = DynRankView("measure", n_cells, n_qp);
            }

            static void weights_times_measure(const DynRankView &q_weights,
                                              const Scalar &factor,
                                              DynRankView &measure) {
                SizeType n_cells = measure.extent(0);
                SizeType n_qp = measure.extent(1);

                assert(n_qp == SizeType(q_weights.extent(0)));

                Kokkos::parallel_for(
                    "ShellTools::weights_times_measure",
                    CellQPRange({0, 0}, {n_cells, n_qp}),
                    KOKKOS_LAMBDA(const SizeType cell, const SizeType qp) {
                        auto &m = measure(cell, qp);
                        auto qw = q_weights(qp);
                        auto w = (factor * qw);
                        m *= w;
                    });
            }

            static void print_ref_grad(const DynRankView &ref_grad) {
                for (int n = 0; n < ref_grad.extent(0); ++n) {
                    for (int qp = 0; qp < ref_grad.extent(1); ++qp) {
                        for (int d_ref = 0; d_ref < ref_grad.extent(2); ++d_ref) {
                            utopia::out() << ref_grad(n, qp, d_ref) << " ";
                        }
                        utopia::out() << "\n";
                    }
                }
            }

            static bool check_quadrature(const DynRankView &q_weights) {
                Scalar wq = sum_of_weights(q_weights);
                utopia::out() << "Weight quadrature: " << wq << "\n";
                return true;
            }

            static Scalar sum_of_weights(const DynRankView &q_weights) {
                typename DynRankView::HostMirror q_weights_host = Kokkos::create_mirror_view(q_weights);
                Scalar wq = 0.0;
                SizeType n_quad_points = q_weights_host.extent(0);
                for (SizeType qp = 0; qp < n_quad_points; ++qp) {
                    wq += q_weights_host(qp);
                }

                return wq;
            }

            // static Scalar ref_measure(const SizeType manifold_dim) {
            //     Scalar ret = 1.0;

            //     // TODO more generic if implemented as 1/(d!)
            //     switch (manifold_dim) {
            //         case 1: {
            //             ret = 1.;
            //             break;
            //         }
            //         case 2: {
            //             ret = 0.5;
            //             break;
            //         }
            //         case 3: {
            //             ret = 1. / 6;
            //             break;
            //         }
            //         default: {
            //             assert(false);
            //         }
            //     }

            //     return ret;
            // }

            // static Scalar ref_grad_scaling(const SizeType manifold_dim) {
            //     Scalar ret = 1.0;

            //     // TODO more generic if implemented as 1/(d!)
            //     switch (manifold_dim) {
            //         case 1: {
            //             ret = 1.;
            //             break;
            //         }
            //         case 2: {
            //             ret = 1.;
            //             break;
            //         }
            //         default: {
            //             ret = 1.0;
            //         }
            //     }

            //     return ret;
            // }

            static void cell_geometry(const DynRankView &cell_nodes,
                                      const DynRankView &q_weights,
                                      const DynRankView &ref_grad,
                                      DynRankView &jacobian,
                                      DynRankView &jacobian_inv,
                                      DynRankView &measure) {
                IsoParametricComputeJacobian<DynRankView> compute_J(cell_nodes, ref_grad, jacobian);
                ComputeJacobianInverse<DynRankView> compute_J_inv(jacobian, jacobian_inv);
                ComputeJacobianDeterminant<DynRankView> compute_J_det(jacobian, measure);

                SizeType n_cells = jacobian.extent(0);
                SizeType n_qp = jacobian.extent(1);

                // Scalar sw = sum_of_weights(q_weights);

                Kokkos::parallel_for(
                    "ShellTools::cell_geometry",
                    CellQPRange({0, 0}, {n_cells, n_qp}),
                    KOKKOS_LAMBDA(const SizeType cell, const SizeType qp) {
                        compute_J(cell, qp);
                        compute_J_det(cell, qp);  // This depends on compute_J(cell)
                        compute_J_inv(cell, qp);  // This depends on compute_J_det(cell)
                    });

                // assert(check_quadrature(q_weights));

                weights_times_measure(q_weights, 1., measure);
            }

            static void transform_gradient_to_physical_space(const DynRankView &jacobian_inv,
                                                             const DynRankView &ref_grad,
                                                             DynRankView &grad_physical) {
                const SizeType n_cells = jacobian_inv.extent(0);
                const SizeType n_qp = jacobian_inv.extent(1);
                const SizeType num_nodes = ref_grad.extent(0);
                const SizeType spatial_dim = grad_physical.extent(3);

                TransformGradientToPhysicalSpace<DynRankView> transform(jacobian_inv, ref_grad, grad_physical);

                // print_ref_grad(ref_grad);

                Kokkos::parallel_for("ShellTools::transform_gradient_to_physical_space",
                                     Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>(
                                         {0, 0, 0, 0}, {n_cells, num_nodes, n_qp, spatial_dim}),
                                     transform);

                // check_grad(grad_physical);
            }
        };  // namespace intrepid2
    }       // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_SHELL_TOOLS_HPP
