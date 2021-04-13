#ifndef UTOPIA_INTREPID2_FORCING_FUNCTION_HPP
#define UTOPIA_INTREPID2_FORCING_FUNCTION_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <typename Fun>
    class ForcingFunction : public Configurable {
    public:
        void read(Input &in) override {
            in.get("value", value);
            in.get("component", component);
        }

        ForcingFunction(const Fun &value = Fun(0.0)) : value(value) {}
        UTOPIA_FUNCTION ForcingFunction(const ForcingFunction &) = default;

        Fun value;
        int n_components{1};
        int component{0};
    };

    namespace intrepid2 {

        template <typename Fun>
        class Assemble<ForcingFunction<Fun>, typename Traits<Fun>::Scalar>
            : public FEAssembler<typename Traits<Fun>::Scalar> {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::ForcingFunction<Fun>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return op_.n_components; }

            inline std::string name() const override { return "ForcingFunction"; }

            bool assemble() override {
                this->ensure_vec_accumulator();

                auto &fe = this->fe();
                auto ev = this->data();

                const int n_components = op_.n_components;
                const int component = op_.component;
                const int n_qp = fe.num_qp();

                auto value = op_.value;
                auto fun = fe.fun;
                auto measure = fe.measure;

                this->vec_integrate(
                    "Assemble<ForcingFunction>::init", KOKKOS_LAMBDA(const int &cell, const int &i) {
                        auto offset = i * n_components + component;

                        for (int qp = 0; qp < n_qp; ++qp) {
                            auto dX = measure(cell, qp);
                            ev(cell, offset) += fun(i, qp) * value * dX;
                        }
                    });

                return true;
            }

            // NVCC_PRIVATE :
            Op op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FORCING_FUNCTION_HPP
