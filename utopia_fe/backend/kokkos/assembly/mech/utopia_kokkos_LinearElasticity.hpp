#ifndef UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP
#define UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"

#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_kokkos_Strain.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

#include "utopia_kokkos_ElasticMaterial.hpp"

#include "utopia_kokkos_PrincipalStresses.hpp"
#include "utopia_kokkos_StressLinearElasticityOp.hpp"

namespace utopia {

    namespace kokkos {

        template <class FE_,
                  int Dim,
                  class FirstLameParameter = typename FE_::Scalar,
                  class ShearModulus = FirstLameParameter>
        class LinearElasticity : public ElasticMaterial<FE_> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::ElasticMaterial<FE_>;
            using VectorView = typename Super::VectorView;

            using Measure = typename FE::Measure;
            using Gradient = typename FE::Gradient;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
                    ssp.read(in);

                    // in.get("mass_density", [&](Input &node) { mass_density.read(node); });

                    lambda = ssp.first_lame_parameter.get();
                    mu = ssp.shear_modulus.get();
                }

                Params(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                       const ShearModulus &mu = FirstLameParameter(1.0))
                    : lambda(lambda), mu(mu) /*, mass_density(1) */ {}

                FirstLameParameter lambda;
                ShearModulus mu;
                // kokkos::SubdomainValue<FE> mass_density;
            };

            using Op = utopia::kokkos::kernels::
                LinearElasticityOp<Dim, Scalar, FirstLameParameter, ShearModulus, Gradient, Measure>;

            LinearElasticity(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline std::string name() const override { return "LinearElasticity"; }

            inline Op make_op() const { return Op(op_.lambda, op_.mu, this->fe().grad(), this->fe().measure()); }

            bool apply(const VectorView &x, VectorView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("LinearElasticity::apply");

                this->apply_vector_operator("LinearElasticity::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("LinearElasticity::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("LinearElasticity::assemble_matrix");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("LinearElasticity::assemble_matrix",
                                           block_op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("LinearElasticity::assemble_matrix");
                return true;
            }

            bool principal_stresses(const VectorView &displacement, VectorView &result) override {
                ////////////////////////////////////////////////////////////////////////////
                using GradientOp = typename utopia::kokkos::kernels::GradientOp<Scalar, Gradient, VectorView>::Rank2;

                using StressOp = utopia::kokkos::kernels::
                    InlineCoeffLinearElasticityOp<Dim, Scalar, FirstLameParameter, ShearModulus, GradientOp, Measure>;

                using PrincipalStressesOp =
                    utopia::kokkos::StorePrincipalStress<Dim, Scalar, StressOp, Measure, VectorView>;
                ////////////////////////////////////////////////////////////////////////////

                UTOPIA_TRACE_REGION_BEGIN("LinearElasticity::principal_stresses");

                if (result.extent(0) != this->fe().n_cells() || int(result.extent(1)) != Dim) {
                    result = VectorView("principal_stresses", this->fe().n_cells(), Dim);
                }

                GradientOp g(this->fe().grad(), displacement);

                StressOp stress(op_.lambda, op_.mu, g, 0);
                PrincipalStressesOp op(stress, this->fe().measure(), result);

                this->loop_cell("LinearElasticity::principal_stresses", op);

                UTOPIA_TRACE_REGION_END("LinearElasticity::principal_stresses");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP
