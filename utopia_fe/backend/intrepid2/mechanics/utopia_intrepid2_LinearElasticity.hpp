#ifndef UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
#define UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_intrepid2_Strain.hpp"
#include "utopia_kokkos_LinearElasticityOp.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class LinearElasticity : public Configurable {
    public:
        void read(Input &in) override {
            StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
            ssp.read(in);

            lambda = ssp.first_lame_parameter.get();
            mu = ssp.shear_modulus.get();
        }

        LinearElasticity(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                         const ShearModulus &mu = FirstLameParameter(1.0))
            : lambda(lambda), mu(mu) {}

        FirstLameParameter lambda;
        ShearModulus mu;
    };

    namespace intrepid2 {

        template <int Dim, class FirstLameParameter, class ShearModulus, typename Scalar>
        class Assemble<LinearElasticity<Dim, FirstLameParameter, ShearModulus>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::LinearElasticity<Dim, FirstLameParameter, ShearModulus>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            using Op = utopia::kokkos::kernels::LinearElasticityOp<Dim,
                                                                   Scalar,
                                                                   FirstLameParameter,
                                                                   ShearModulus,
                                                                   typename FE::Gradient,
                                                                   typename FE::Measure>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline std::string name() const override { return "LinearElasticity"; }

            inline Op make_op() const { return Op(op_.lambda, op_.mu, this->fe().grad(), this->fe().measure()); }

            bool apply(const DynRankView &x, DynRankView &y) override {
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
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
