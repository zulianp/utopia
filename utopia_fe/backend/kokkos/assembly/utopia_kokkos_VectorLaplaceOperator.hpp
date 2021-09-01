#ifndef UTOPIA_KOKKOS_VECTOR_LAPLACE_OPERATOR_HPP
#define UTOPIA_KOKKOS_VECTOR_LAPLACE_OPERATOR_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_LaplaceOp.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, int Dim, typename Coefficient = typename FE_::Scalar>
        class VectorLaplaceOperator : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("coeff", coeff);
                    in.get("verbose", verbose);
                }

                Params(const Coefficient &coeff = Coefficient(1.0)) : coeff(coeff) {}
                Coefficient coeff;
                bool verbose{false};
            };

            using Op = utopia::kokkos::kernels::
                VectorLaplaceOp<Dim, Scalar, Coefficient, typename FE::Gradient, typename FE::Measure>;

            VectorLaplaceOperator(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }
            inline std::string name() const override { return "VectorLaplaceOperator"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline Op make_op() const { return Op(op_.coeff, this->fe().grad(), this->fe().measure()); }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("VectorLaplaceOperator::apply");

                this->apply_vector_operator("VectorLaplaceOperator::apply", x, y, make_op());

                if (op_.verbose) {
                    utopia::out() << "ForcingFunction: " << this->vector_accumulator()->sum() << '\n';
                    // this->describe(utopia::out().stream());
                    // utopia::out() << "Accumulator:\n";
                    this->vector_accumulator()->describe(utopia::out().stream());
                }

                UTOPIA_TRACE_REGION_END("VectorLaplaceOperator::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("VectorLaplaceOperator::assemble_matrix");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("VectorLaplaceOperator::assemble_matrix",
                                           block_op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("VectorLaplaceOperator::assemble_matrix");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_VECTOR_LAPLACE_OPERATOR_HPP
