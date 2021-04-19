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
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::Transport<Dim, Field>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : Super(fe), op_(std::move(op)) {
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

            int rank() const override { return 2; }
            inline std::string name() const override { return "Transport"; }
            bool assemble() override {
                this->ensure_mat_accumulator();

                auto &fe = this->fe();
                auto data = this->data();

                const int n_qp = fe.num_qp();

                {
                    auto vector_field = op_.vector_field;
                    auto grad = fe.grad;
                    auto fun = fe.fun;
                    auto measure = fe.measure;

                    this->mat_integrate(
                        "Assemble<Transport>::init", KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            assert((Dim == 1 || &grad(cell, j, 0, 1) - &grad(cell, j, 0, 0) == 1UL) &&
                                   "spatial dimension must be contiguous");

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

                            data(cell, i, j) += integral;
                        });
                }

                return true;
            }

            // NVCC_PRIVATE :
            Op op_;
        };
    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_TRANSPORT_HPP