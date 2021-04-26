#ifndef UTOPIA_INTREPID2_NEOHOOKEAN_HPP
#define UTOPIA_INTREPID2_NEOHOOKEAN_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Gradient.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class NeoHookean : public Configurable {
    public:
        using Scalar = typename Traits<FirstLameParameter>::Scalar;
        void read(Input &in) override {
            in.get("lambda", lambda);
            in.get("mu", mu);
        }

        NeoHookean(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                   const ShearModulus &mu = FirstLameParameter(1.0),
                   const Scalar &rescale = 1.0)
            : lambda(lambda), mu(mu), rescale(rescale) {}

        FirstLameParameter lambda;
        ShearModulus mu;
        Scalar rescale;
    };

    namespace intrepid2 {

        template <int Dim, class FirstLameParameter, class ShearModulus, typename Scalar>
        class Assemble<NeoHookean<Dim, FirstLameParameter, ShearModulus>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::NeoHookean<Dim, FirstLameParameter, ShearModulus>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }

            inline std::string name() const override { return "NeoHookean"; }

            virtual bool update(const std::shared_ptr<Field<Scalar>> &displacement) {
                assert(displacement);
                assert(displacement->is_coefficient());

                if (!deformation_gradient_) {
                    // Initialize gradient
                    deformation_gradient_ = std::make_shared<Gradient<Scalar>>(this->fe());
                }

                deformation_gradient_->init(*displacement);
            }

            class OpAndStoreMat {
            public:
                using V = StaticVector<Scalar, Dim>;
                using M = StaticMatrix<Scalar, Dim, Dim>;

                OpAndStoreMat(const FirstLameParameter &lambda,
                              const ShearModulus &mu,
                              const Scalar &rescaling,
                              const DynRankView &F,
                              const DynRankView &grad,
                              const DynRankView &measure,
                              DynRankView &data)
                    : lambda(lambda), mu(mu), rescaling(rescaling), F(F), grad(grad), measure(measure), data(data) {}

                inline static void selective_triple_product(const int di,
                                                            const M &F_inv_t_qp,
                                                            const V &grad,
                                                            M &FGF_i) {
                    M grad_mat;
                    // TODO Make more efficient by skipping 0 operations
                    make_grad(di, grad, grad_mat);
                    FGF_i = F_inv_t_qp * transpose(grad_mat) * F_inv_t_qp;
                }

                inline static Scalar selective_inner(const int di, const M &mat, const V &grad) {
                    Scalar ret = 0.0;

                    for (int dj = 0; dj < Dim; ++dj) {
                        ret += mat(di, dj) * grad[dj];
                    }

                    return ret;
                }

                inline static void make_grad(const int di, const V &vec_grad, M &mat_grad) {
                    mat_grad.set(0.0);

                    for (int i = 0; i < Dim; ++i) {
                        mat_grad(di, i) = vec_grad[i];
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const int &cell, const int qp, const int &i, const int &j) {
                    V grad_i, grad_j;
                    M FGF_i, stress_lin_i;
                    M F_qp, F_inv_t_qp;

                    auto dX = measure(cell, qp);

                    for (int d = 0; d < Dim; ++d) {
                        grad_i[d] = grad(cell, i, qp, d);
                        grad_j[d] = grad(cell, j, qp, d);
                    }

                    for (int d1 = 0; d1 < Dim; ++d1) {
                        for (int d2 = 0; d2 < Dim; ++d2) {
                            F_qp(d1, d2) = F(cell, qp, d1, d2);
                        }
                    }

                    F_inv_t_qp = inv(transpose(F_qp));
                    Scalar J = det(F_qp);
                    Scalar log_J = device::log(J);
                    Scalar beta = (rescaling * (lambda * log_J - mu));

                    for (int di = 0; di < Dim; ++di) {
                        auto dof_i = i * Dim + di;

                        selective_triple_product(di, F_inv_t_qp, grad_i, FGF_i);
                        Scalar inner_F_grad_i = selective_inner(di, F_inv_t_qp, grad_i);

                        make_grad(di, grad_i, stress_lin_i);
                        stress_lin_i *= rescaling * mu;
                        stress_lin_i += beta * FGF_i + (lambda * inner_F_grad_i) * F_inv_t_qp;

                        for (int dj = 0; dj < Dim; ++dj) {
                            auto dof_j = j * Dim + dj;
                            auto val = selective_inner(dj, stress_lin_i, grad_j);
                            data(cell, dof_i, dof_j) += val;
                        }
                    }
                }

                const FirstLameParameter lambda;
                const ShearModulus mu;
                const Scalar rescaling;
                const DynRankView F;
                const DynRankView grad;
                const DynRankView measure;
                DynRankView data;
            };

            bool assemble_jacobian() {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::assemble_jacobian");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    auto F = deformation_gradient_->data();
                    this->loop_cell_test_trial(
                        "Assemble<NeoHookean>::assemble_jacobian",
                        OpAndStoreMat(op_.lambda, op_.mu, op_.rescaling, F, fe.grad, fe.measure, data));
                }

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::assemble_jacobian");
            }

            bool assemble() override { return assemble_jacobian(); }

            // NVCC_PRIVATE :
            UserOp op_;
            // std::shared_ptr<Field<Scalar>> displacement_;
            std::shared_ptr<Gradient<Scalar>> deformation_gradient_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_NEOHOOKEAN_HPP
