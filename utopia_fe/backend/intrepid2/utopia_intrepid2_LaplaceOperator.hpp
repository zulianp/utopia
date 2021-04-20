#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

namespace utopia {

    template <class DiffusionCoefficient>
    class LaplaceOperator : public Configurable {
    public:
        void read(Input &in) override { in.get("coeff", coeff); }

        LaplaceOperator(const DiffusionCoefficient &coeff = DiffusionCoefficient(1.0)) : coeff(coeff) {}
        DiffusionCoefficient coeff;
    };

    namespace intrepid2 {

        template <class DiffusionCoefficient, typename Scalar>
        class Assemble<LaplaceOperator<DiffusionCoefficient>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::LaplaceOperator<DiffusionCoefficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            int rank() const override { return 2; }
            inline std::string name() const override { return "LaplaceOperator"; }

            bool assemble() override {
                this->ensure_mat_accumulator();

                auto &fe = this->fe();
                auto data = this->data();

                DynRankView grad_x_measure(
                    "grad_x_measure", fe.num_cells(), fe.num_fields(), fe.num_qp(), fe.spatial_dimension());

                Kokkos::deep_copy(grad_x_measure, fe.grad);
                FunctionSpaceTools::template multiplyMeasure<Scalar>(grad_x_measure, fe.measure, fe.grad);
                FunctionSpaceTools::template integrate<Scalar>(data, fe.grad, grad_x_measure);

                // Only works if coeff is a scalar
                {
                    auto c = op_.coeff;
                    this->loop_cell_test_trial(
                        "scale_with_coeff",
                        KOKKOS_LAMBDA(const int &i0, const int &i1, const int &i2) { data(i0, i1, i2) *= c; });
                }

                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif
