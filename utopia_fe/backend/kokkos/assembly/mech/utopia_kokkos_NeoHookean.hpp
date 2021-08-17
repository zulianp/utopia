#ifndef UTOPIA_KOKKOS_NEOHOOKEAN_HPP
#define UTOPIA_KOKKOS_NEOHOOKEAN_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_NeoHookeanOp.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, class FirstLameParameter = typename FE_::Scalar, class ShearModulus = FirstLameParameter>
        class NeoHookean : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;

            template <int Dim>
            struct Def {
                using NeoHookeanKernel = utopia::kokkos::kernels::NeoHookeanOp<Dim,
                                                                               Scalar,
                                                                               FirstLameParameter,
                                                                               ShearModulus,
                                                                               DynRankView,
                                                                               typename FE::Gradient,
                                                                               typename FE::Measure>;

                using Op = typename NeoHookeanKernel::ApplyHessian;
                using OpAndStoreHessian = typename NeoHookeanKernel::template StoreHessian<DynRankView>;
                using OpAndStoreGradient = typename NeoHookeanKernel::template StoreGradient<DynRankView>;
            };

            class Params : public Configurable {
            public:
                using Scalar = typename Traits<FirstLameParameter>::Scalar;
                void read(Input &in) override {
                    StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
                    ssp.read(in);

                    lambda = ssp.first_lame_parameter.get();
                    mu = ssp.shear_modulus.get();
                }

                Params(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                       const ShearModulus &mu = FirstLameParameter(1.0),
                       const Scalar &rescale = 1.0)
                    : lambda(lambda), mu(mu), rescale(rescale) {}

                FirstLameParameter lambda;
                ShearModulus mu;
                Scalar rescale;
            };

            NeoHookean(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return this->fe().spatial_dimension(); }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "NeoHookean"; }

            virtual bool update(const std::shared_ptr<Field<FE>> &displacement) override {
                if (!Super::update(displacement)) {
                    return false;
                }

                assert(displacement);
                assert(displacement->is_coefficient());

                if (!displacement->is_coefficient()) {
                    Utopia::Abort("Assemble<NeoHookean>::update, displacement must me in coefficient form!");
                }

                if (!deformation_gradient_) {
                    // Initialize gradient
                    deformation_gradient_ = std::make_shared<Gradient<FE>>(this->fe_ptr());
                }

                deformation_gradient_->init(*displacement);
                deformation_gradient_->add_identity();
                assert(deformation_gradient_->check_dets_are_positive());
                return true;
            }

            template <int Dim>
            inline typename Def<Dim>::Op make_op() const {
                assert(deformation_gradient_);
                auto F = deformation_gradient_->data();
                return
                    typename Def<Dim>::Op(op_.lambda, op_.mu, op_.rescale, F, this->fe().grad(), this->fe().measure());
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::apply");

                switch (this->fe().spatial_dimension()) {
                    case 1: {
                        this->apply_vector_operator("Assemble<NeoHookean>::apply", x, y, make_op<1>());
                        break;
                    }
                    case 2: {
                        this->apply_vector_operator("Assemble<NeoHookean>::apply", x, y, make_op<2>());
                        break;
                    }
                    case 3: {
                        this->apply_vector_operator("Assemble<NeoHookean>::apply", x, y, make_op<3>());
                        break;
                    }
                    default: {
                        assert(false);
                        Utopia::Abort("Invalid dimension!");
                    }
                }

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();

                    switch (this->fe().spatial_dimension()) {
                        case 1: {
                            this->loop_cell_test_trial(
                                "Assemble<NeoHookean>::assemble_matrix",
                                typename Def<1>::OpAndStoreHessian(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));
                            break;
                        }

                        case 2: {
                            this->loop_cell_test_trial(
                                "Assemble<NeoHookean>::assemble_matrix",
                                typename Def<2>::OpAndStoreHessian(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));
                            break;
                        }

                        case 3: {
                            this->loop_cell_test_trial(
                                "Assemble<NeoHookean>::assemble_matrix",
                                typename Def<3>::OpAndStoreHessian(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));
                            break;
                        }

                        default: {
                            assert(false);
                            Utopia::Abort("Invalid dimension!");
                        }
                    }
                }

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::assemble_vector");

                this->ensure_vector_accumulator();

                auto &fe = this->fe();
                auto data = this->vector_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();
                    // this->loop_cell_test(
                    //     "Assemble<NeoHookean>::assemble_vector",
                    //     OpAndStoreGradient(op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));

                    switch (this->fe().spatial_dimension()) {
                        case 1: {
                            this->loop_cell_test(
                                "Assemble<NeoHookean>::assemble_vector",
                                typename Def<1>::OpAndStoreGradient(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));

                            break;
                        }

                        case 2: {
                            this->loop_cell_test(
                                "Assemble<NeoHookean>::assemble_vector",
                                typename Def<3>::OpAndStoreGradient(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));

                            break;
                        }

                        case 3: {
                            this->loop_cell_test(
                                "Assemble<NeoHookean>::assemble_vector",
                                typename Def<3>::OpAndStoreGradient(
                                    op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));

                            break;
                        }

                        default: {
                            assert(false);
                            Utopia::Abort("Invalid dimension!");
                        }
                    }
                }

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
            std::shared_ptr<Gradient<FE>> deformation_gradient_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_NEOHOOKEAN_HPP
