#ifndef UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HPP
#define UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HPP

#include "utopia_kokkos_AutoHyperElasticity.hpp"
#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

#include "generated_AutoHyperElasticity.hpp"
#include "generated_AutoHyperElasticityEnergy.hpp"
#include "generated_AutoHyperElasticityGradient.hpp"
#include "generated_AutoHyperElasticityHessian.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, class FirstLameParameter = typename FE_::Scalar, class ShearModulus = FirstLameParameter>
        class AutoHyperElasticity : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;

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

            AutoHyperElasticity(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return this->fe().spatial_dimension(); }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "AutoHyperElasticity"; }

            virtual bool update(const std::shared_ptr<Field<FE>> &displacement) override {
                if (!Super::update(displacement)) {
                    return false;
                }

                assert(displacement);
                assert(displacement->is_coefficient());

                if (!displacement->is_coefficient()) {
                    Utopia::Abort("Assemble<AutoHyperElasticity>::update, displacement must me in coefficient form!");
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

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoHyperElasticity>::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();
                }

                Utopia::Abort("IMPLEMENT ME");

                UTOPIA_TRACE_REGION_END("Assemble<AutoHyperElasticity>::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<AutoHyperElasticity>::assemble_vector");

                this->ensure_vector_accumulator();

                auto &fe = this->fe();
                auto data = this->vector_data();

                constexpr int Dim = 2;

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();

                    auto lambda = op_.lambda;
                    auto mu = op_.mu;

                    int n_quad_points = this->fe().n_quad_points();
                    int n_shape_functions = this->fe().n_shape_functions();

                    auto grad = this->fe().grad();
                    auto measure = this->fe().measure();

                    this->loop_cell(
                        "Gradient", UTOPIA_LAMBDA(int cell) {
                            StaticVector<Scalar, Dim> stress, grad_test;
                            StaticMatrix<Scalar, Dim, Dim> F_qp;

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                // Copy deformation gradient at qp
                                for (int d1 = 0; d1 < Dim; ++d1) {
                                    for (int d2 = 0; d2 < Dim; ++d2) {
                                        F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                                    }
                                }

                                Scalar dx = measure(cell, qp);

                                for (int i = 0; i < n_shape_functions; ++i) {
                                    // Copy deformation gradient at qp
                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        grad_test[d1] = grad(cell, i, qp, d1);
                                    }

                                    // Assemble
                                    stress.set(0.);
                                    elastic_material_gradient(mu, lambda, &F_qp.raw_type()[0], &grad_test[0], dx, &stress.raw_type()[0], 0);

                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        auto dof_i = i * Dim + d1;
                                        data(cell, dof_i) += stress[d1];
                                    }
                                }
                            }
                        });
                }

                UTOPIA_TRACE_REGION_END("Assemble<AutoHyperElasticity>::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
            std::shared_ptr<Gradient<FE>> deformation_gradient_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HPP
