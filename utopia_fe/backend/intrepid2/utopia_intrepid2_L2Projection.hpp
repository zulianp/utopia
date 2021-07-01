#ifndef UTOPIA_INTREPID2_L2_PROJECTION_HPP
#define UTOPIA_INTREPID2_L2_PROJECTION_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Gradient.hpp"

#include "utopia_Views.hpp"

#include <string>

namespace utopia {

    template <class Field>
    class L2Projection : public Configurable {
    public:
        void read(Input &) override {}

        L2Projection() = default;
        L2Projection(const Field &field) : field(field) {}
        Field field;
    };

    namespace intrepid2 {

        template <class Field>
        class Assemble<L2Projection<Field>, typename Traits<Field>::Scalar>
            : public utopia::intrepid2::FEAssembler<typename Traits<Field>::Scalar> {
        public:
            using Scalar = typename Traits<Field>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using UserOp = utopia::L2Projection<Field>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<typename Traits<Field>::Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, const UserOp &op) : Super(fe), op_(op) {}

            inline int n_vars() const override { return op_.field.extent(2); }

            inline std::string name() const override { return "L2Projection"; }

            inline bool is_matrix() const override { return false; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return false; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Field &field, const DynRankView &fun, const DynRankView &measure)
                    : field(field), fun(fun), measure(measure), n_qp(measure.extent(1)) {
                    assert(field.extent(0) == measure.extent(0));
                    assert(field.extent(1) == measure.extent(1));
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int sub_i) const {
                    Scalar ret = 0;

                    for (int qp = 0; qp < n_qp; ++qp) {
                        ret += field(cell, qp, sub_i) * fun(i, qp) * measure(qp);
                    }

                    return ret;
                }

                Field field;
                DynRankView fun;
                DynRankView measure;
                const int n_qp;
            };

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.field, fe.fun, fe.measure);
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<L2Projection>::assemble");

                this->ensure_vector_accumulator();
                this->loop_cell_test_trial("Assemble<L2Projection>::assemble",
                                           op_and_store_cell_ij(this->vector_data(), make_op()));

                UTOPIA_TRACE_REGION_END("Assemble<L2Projection>::assemble");
                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_L2_PROJECTION_HPP
