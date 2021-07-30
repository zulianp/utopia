#ifndef UTOPIA_INTREPID2_NEOHOOKEAN_HPP
#define UTOPIA_INTREPID2_NEOHOOKEAN_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Gradient.hpp"
#include "utopia_kokkos_NeoHookeanOp.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class NeoHookean : public Configurable {
    public:
        using Scalar = typename Traits<FirstLameParameter>::Scalar;
        void read(Input &in) override {
            StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
            ssp.read(in);

            lambda = ssp.first_lame_parameter.get();
            mu = ssp.shear_modulus.get();
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

            using V = StaticVector<Scalar, Dim>;
            using M = StaticMatrix<Scalar, Dim, Dim>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "NeoHookean"; }

            virtual bool update(const std::shared_ptr<Field<Scalar>> &displacement) override {
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
                    deformation_gradient_ = std::make_shared<Gradient<Scalar>>(this->fe_ptr());
                }

                deformation_gradient_->init(*displacement);
                deformation_gradient_->add_identity();
                assert(deformation_gradient_->check_dets_are_positive());
                return true;
            }

            inline Op make_op() const {
                assert(deformation_gradient_);
                auto F = deformation_gradient_->data();

                return Op(op_.lambda, op_.mu, op_.rescale, F, this->fe().grad(), this->fe().measure());
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::apply");

                this->apply_vector_operator("Assemble<NeoHookean>::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::apply");
                return false;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<NeoHookean>::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();
                    this->loop_cell_test_trial(
                        "Assemble<NeoHookean>::assemble_matrix",
                        OpAndStoreHessian(op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));
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
                    this->loop_cell_test(
                        "Assemble<NeoHookean>::assemble_vector",
                        OpAndStoreGradient(op_.lambda, op_.mu, op_.rescale, F, fe.grad(), fe.measure(), data));
                }

                UTOPIA_TRACE_REGION_END("Assemble<NeoHookean>::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
            std::shared_ptr<Gradient<Scalar>> deformation_gradient_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_NEOHOOKEAN_HPP
