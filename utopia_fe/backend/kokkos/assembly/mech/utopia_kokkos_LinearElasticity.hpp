#ifndef UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP
#define UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"

#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_kokkos_Strain.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    namespace kokkos {

        template <class FE_,
                  int Dim,
                  class FirstLameParameter = typename FE_::Scalar,
                  class ShearModulus = FirstLameParameter>
        class LinearElasticity : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;
            using VectorView = typename Super::VectorView;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
                    ssp.read(in);

                    lambda = ssp.first_lame_parameter.get();
                    mu = ssp.shear_modulus.get();
                }

                Params(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                       const ShearModulus &mu = FirstLameParameter(1.0))
                    : lambda(lambda), mu(mu) {}

                FirstLameParameter lambda;
                ShearModulus mu;
            };

            using Op = utopia::kokkos::kernels::LinearElasticityOp<Dim,
                                                                   Scalar,
                                                                   FirstLameParameter,
                                                                   ShearModulus,
                                                                   typename FE::Gradient,
                                                                   typename FE::Measure>;

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
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LinearElasticity>::apply");

                this->apply_vector_operator("Assemble<LinearElasticity>::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("Assemble<LinearElasticity>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LinearElasticity>::assemble_matrix");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("Assemble<LinearElasticity>::assemble_matrix",
                                           block_op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("Assemble<LinearElasticity>::assemble_matrix");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LINEAR_ELASTICITY_HPP
