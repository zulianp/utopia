#ifndef UTOPIA_INTREPID_2_MASS_HPP
#define UTOPIA_INTREPID_2_MASS_HPP

#include "utopia_intrepid2_LaplaceOperator.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <typename Fun>
    class Mass : public Configurable {
    public:
        void read(Input &in) override {
            in.get("density", density);
            in.get("n_components", n_components);
        }

        Mass(const Fun &density = Fun(1.0)) : density(density) {}
        UTOPIA_FUNCTION Mass(const Mass &) = default;

        Fun density;
        int n_components{1};
    };

    namespace intrepid2 {

        template <typename Fun>
        class Assemble<Mass<Fun>, typename Traits<Fun>::Scalar> : public FEAssembler<typename Traits<Fun>::Scalar> {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::Mass<Fun>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return op_.n_components; }
            inline std::string name() const override { return "Mass"; }

            bool assemble() override {
                this->ensure_mat_accumulator();
                auto &fe = this->fe();
                auto data = this->data();

                const int n_components = op_.n_components;
                // const int component = op_.component;
                const int n_qp = fe.num_qp();

                {
                    auto density = op_.density;
                    auto fun = fe.fun;
                    auto measure = fe.measure;

                    assert(n_components == 1 && "IMPLEMENT ME");

                    this->mat_integrate(
                        "Assemble<Mass>::init", KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            auto offset_i = i * n_components;
                            auto offset_j = j * n_components;

                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);
                                data(cell, offset_i, offset_j) += fun(i, qp) * fun(j, qp) * density * dX;
                            }
                        });
                }

                return true;
            }

            // NVCC_PRIVATE :
            Op op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID_2_MASS_HPP
