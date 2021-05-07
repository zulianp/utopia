#ifndef UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class Coefficient>
    class VectorLaplaceOperator : public Configurable {
    public:
        void read(Input &in) override { in.get("coeff", coeff); }

        VectorLaplaceOperator(const Coefficient &coeff = Coefficient(1.0)) : coeff(coeff) {}
        Coefficient coeff;
    };

    namespace intrepid2 {

        template <int Dim, typename Coefficient, typename Scalar>
        class Assemble<VectorLaplaceOperator<Dim, Coefficient>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::VectorLaplaceOperator<Dim, Coefficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }
            inline std::string name() const override { return "VectorLaplaceOperator"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return false; }
            inline bool is_scalar() const override { return false; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;

                        for (int d = 0; d < dim(); ++d) {
                            val += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += coeff * val * measure(cell, qp);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                UTOPIA_INLINE_FUNCTION Op(const Coefficient &coeff, const DynRankView &grad, const DynRankView &measure)
                    : coeff(coeff), grad(grad), measure(measure), n_qp(measure.extent(1)) {}

                const Coefficient coeff;
                const DynRankView grad;
                const DynRankView measure;
                const int n_qp;
            };

            inline Op make_op() const { return Op(op_.coeff, this->fe().grad, this->fe().measure); }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<VectorLaplaceOperator>::apply");

                this->apply_vector_operator("Assemble<VectorLaplaceOperator>::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("Assemble<VectorLaplaceOperator>::apply");
                return false;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<VectorLaplaceOperator>::assemble_matrix");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("Assemble<VectorLaplaceOperator>::assemble_matrix",
                                           block_op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("Assemble<VectorLaplaceOperator>::assemble_matrix");
                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
