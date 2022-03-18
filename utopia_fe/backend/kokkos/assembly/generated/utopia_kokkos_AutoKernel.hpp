#ifndef UTOPIA_ASSEMBLY_KOKKOS_AUTOKERNEL_HPP
#define UTOPIA_ASSEMBLY_KOKKOS_AUTOKERNEL_HPP

#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"

#include "utopia_kokkos_FEAssembler.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, class Kernel, int Dim>
        class AutoKernel {};

        template <class FE_, class Kernel>
        class AutoKernel<FE_, Kernel, 2> : public FEAssembler<FE_> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using Params = typename Kernel::Params;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_>;

            class MatrixKernel {
            public:
                static constexpr int NNodes = Kernel::ElemT::NNodes;

                MatrixKernel(const FE &fe, DynRankView &result)
                    : points(fe.points()), result(result), q_points(fe.q_points), q_weights(fe.q_weights) {}

                UTOPIA_FUNCTION void operator()(const int cell) const {
                    Scalar x[NNodes], y[NNodes];
                    Scalar H[NNodes * NNodes];

                    // Initialize values to zero
                    for (int i = 0; i < NNodes * NNodes; ++i) {
                        H[i] = 0;
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        x[i] = points(cell, i, 0);
                        y[i] = points(cell, i, 1);
                    }

                    const int n_quad_points = q_weights.extent(0);

                    for (int q = 0; q < n_quad_points; ++q) {
                        kernel.hessian(x, y, nullptr /*FIXME*/, q_points(q, 0), q_points(q, 1), q_weights(q), H);
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        for (int j = 0; j < NNodes * NNodes; ++j) {
                            result(cell, i, j) += H[i * NNodes + j];
                        }
                    }
                }

                DynRankView points;
                DynRankView result;

                DynRankView q_points;
                DynRankView q_weights;

                Kernel kernel;
            };

            class ApplyKernel {
            public:
                static constexpr int NNodes = Kernel::ElemT::NNodes;

                ApplyKernel(const FE &fe, const DynRankView &in, DynRankView &result)
                    : points(fe.points()), in(in), result(result), q_points(fe.q_points), q_weights(fe.q_weights) {}

                UTOPIA_FUNCTION void operator()(const int cell) const {
                    Scalar x[NNodes], y[NNodes];
                    Scalar Hx[NNodes];
                    Scalar u[NNodes];

                    // Initialize values to zero
                    for (int i = 0; i < NNodes; ++i) {
                        Hx[i] = 0;
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        u[i] = in(cell, i);
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        x[i] = points(cell, i, 0);
                        y[i] = points(cell, i, 1);
                    }

                    const int n_quad_points = q_weights.extent(0);

                    for (int q = 0; q < n_quad_points; ++q) {
                        kernel.apply(x, y, u, q_points(q, 0), q_points(q, 1), q_weights(q), Hx);
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        result(cell, i) += Hx[i];
                    }
                }

                DynRankView points;
                DynRankView in;
                DynRankView result;

                DynRankView q_points;
                DynRankView q_weights;

                Kernel kernel;
            };

            AutoKernel(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "AutoKernel"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoKernel>::apply");

                {
                    ApplyKernel kernel(this->fe(), x, y);
                    this->loop_cell("ApplyHessian", kernel);
                }

                UTOPIA_TRACE_REGION_END("Assemble<AutoKernel>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoKernel>::assemble");

                this->ensure_matrix_accumulator();

                {
                    auto data = this->matrix_data();
                    MatrixKernel kernel(this->fe(), data);
                    this->loop_cell("Hessian", kernel);
                }

                UTOPIA_TRACE_REGION_END("Assemble<AutoKernel>::assemble");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };

        template <class FE_, class Kernel>
        class AutoKernel<FE_, Kernel, 3> : public FEAssembler<FE_> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using Params = typename Kernel::Params;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_>;

            class MatrixKernel {
            public:
                static constexpr int NNodes = Kernel::ElemT::NNodes;

                MatrixKernel(const FE &fe, DynRankView &result)
                    : points(fe.points()), result(result), q_points(fe.q_points), q_weights(fe.q_weights) {}

                UTOPIA_FUNCTION void operator()(const int cell) const {
                    Scalar x[NNodes], y[NNodes], z[NNodes];
                    Scalar H[NNodes * NNodes];

                    // Initialize values to zero
                    for (int i = 0; i < NNodes * NNodes; ++i) {
                        H[i] = 0;
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        x[i] = points(cell, i, 0);
                        y[i] = points(cell, i, 1);
                        z[i] = points(cell, i, 2);
                    }

                    const int n_quad_points = q_weights.extent(0);

                    Scalar sum_w = 0;

                    for (int q = 0; q < n_quad_points; ++q) {
                        Scalar qx = q_points(q, 0);
                        Scalar qy = q_points(q, 1);
                        Scalar qz = q_points(q, 2);
                        Scalar qw = q_weights(q);

                        sum_w += qw;

                        assert(qx >= 0);
                        assert(qx < 1);

                        kernel.hessian(x, y, z, nullptr /*FIXME*/, qx, qy, qz, qw, H);
                    }

                    Scalar vol = 0;
                    for (int i = 0; i < NNodes; ++i) {
                        for (int j = 0; j < NNodes; ++j) {
                            result(cell, i, j) += H[i * NNodes + j];
                            vol += H[i * NNodes + j];
                        }
                    }
                }

                DynRankView points;
                DynRankView result;

                DynRankView q_points;
                DynRankView q_weights;

                Kernel kernel;
            };

            class ApplyKernel {
            public:
                static constexpr int NNodes = Kernel::ElemT::NNodes;

                ApplyKernel(const FE &fe, const DynRankView &in, DynRankView &result)
                    : points(fe.points()), in(in), result(result), q_points(fe.q_points), q_weights(fe.q_weights) {}

                UTOPIA_FUNCTION void operator()(const int cell) const {
                    Scalar x[NNodes], y[NNodes], z[NNodes];
                    Scalar Hx[NNodes];
                    Scalar u[NNodes];

                    // Initialize values to zero
                    for (int i = 0; i < NNodes; ++i) {
                        Hx[i] = 0;
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        u[i] = in(cell, i);
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        x[i] = points(cell, i, 0);
                        y[i] = points(cell, i, 1);
                        z[i] = points(cell, i, 2);
                    }

                    const int n_quad_points = q_weights.extent(0);

                    for (int q = 0; q < n_quad_points; ++q) {
                        kernel.apply(x, y, z, u, q_points(q, 0), q_points(q, 1), q_points(q, 2), q_weights(q), Hx);
                    }

                    for (int i = 0; i < NNodes; ++i) {
                        result(cell, i) += Hx[i];
                    }
                }

                DynRankView points;
                DynRankView in;
                DynRankView result;

                DynRankView q_points;
                DynRankView q_weights;

                Kernel kernel;
            };

            AutoKernel(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "AutoKernel"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoKernel>::apply");

                {
                    ApplyKernel kernel(this->fe(), x, y);
                    this->loop_cell("ApplyHessian", kernel);
                }

                UTOPIA_TRACE_REGION_END("Assemble<AutoKernel>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoKernel>::assemble");

                this->ensure_matrix_accumulator();

                {
                    auto data = this->matrix_data();
                    MatrixKernel kernel(this->fe(), data);
                    this->loop_cell("Hessian", kernel);
                }

                UTOPIA_TRACE_REGION_END("Assemble<AutoKernel>::assemble");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_ASSEMBLY_KOKKOS_AUTOKERNEL_HPP
