#ifndef UTOPIA_INTREPID2_TRANSPORT_HPP
#define UTOPIA_INTREPID2_TRANSPORT_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Gradient.hpp"

#include "utopia_Views.hpp"

namespace utopia {

    template <int Dim, class Field>
    class Transport : public Configurable {
    public:
        void read(Input &in) override {}

        Transport() = default;
        Transport(const Field &vector_field) : vector_field(vector_field) {}
        Field vector_field;
    };

    namespace intrepid2 {

        template <int Dim, class Field>
        class Assemble<Transport<Dim, Field>, typename Traits<Field>::Scalar>
            : public utopia::intrepid2::FEAssembler<typename Traits<Field>::Scalar> {
        public:
            using Scalar = typename Traits<Field>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using UserOp = utopia::Transport<Dim, Field>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());

#ifndef NDEBUG
                if (op.vector_field.size() > 0) {
                    assert(op.vector_field.extent(0) == fe->num_cells());
                    assert(op.vector_field.extent(1) == fe->num_qp());
                    assert(op.vector_field.extent(2) == fe->spatial_dimension());
                }
#endif  // NDEBUG
            }

            inline int n_vars() const override { return Dim; }

            inline std::string name() const override { return "Transport"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Field &vector_field,
                                          const DynRankView &grad,
                                          const DynRankView &fun,
                                          const DynRankView &measure)
                    : vector_field(vector_field), grad(grad), fun(fun), measure(measure), n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar integral = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);

                        Scalar val = 0.0;
                        for (int dj = 0; dj < Dim; ++dj) {
                            val += grad(cell, j, qp, 0, dj) * vector_field(cell, qp, dj) * dX;
                        }

                        val *= fun(cell, i, qp) * dX;
                        integral += val;
                    }

                    return integral;
                }

                Field vector_field;
                DynRankView grad;
                DynRankView fun;
                DynRankView measure;
                const int n_qp;
            };

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.vector_field, fe.grad, fe.fun, fe.measure);
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<Transport>::apply");

                this->apply_operator("Assemble<Transport>::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("Assemble<Transport>::apply");
                return true;
            }

            bool assemble() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<Transport>::assemble");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("Assemble<Transport>::assemble",
                                           op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("Assemble<Transport>::assemble");
                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_TRANSPORT_HPP
