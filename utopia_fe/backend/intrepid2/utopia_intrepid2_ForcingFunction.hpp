#ifndef UTOPIA_INTREPID2_FORCING_FUNCTION_HPP
#define UTOPIA_INTREPID2_FORCING_FUNCTION_HPP

#include "utopia_intrepid2_LaplaceOperator.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <typename Fun>
    class ForcingFunction : public Configurable {
    public:
        void read(Input &in) override {
            in.get("value", value);
            in.get("component", component);
        }

        ForcingFunction(const Fun &value) : value(value) {}

        UTOPIA_FUNCTION ForcingFunction() = default;
        UTOPIA_FUNCTION ForcingFunction(const ForcingFunction &) = default;

        Fun value;
        int n_components{1};
        int component{0};
    };

    namespace intrepid2 {

        template <typename Fun>
        class Assemble<ForcingFunction<Fun>, typename Traits<Fun>::Scalar> : public Describable {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::ForcingFunction<Fun>;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : op_(std::move(op)), fe_(fe) {}

            void init() {
                const int n_components = op_.n_components;
                const int component = op_.component;

                const int num_fields = fe_->num_fields();
                const int n_dofs = num_fields * n_components;
                const int n_qp = fe_->num_qp();

                element_vectors_ = DynRankView("ForcingFunction", fe_->num_cells(), n_dofs);

                {
                    auto ev = element_vectors_;
                    auto value = op_.value;
                    auto fun = fe_->fun;
                    auto measure = fe_->measure;

                    Kokkos::parallel_for(
                        "Assemble<ForcingFunction>::init",
                        Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>(
                            {0, 0}, {static_cast<int>(ev.extent(0)), num_fields}),
                        KOKKOS_LAMBDA(const int &cell, const int &i) {
                            auto offset = i * n_components + component;

                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);
                                ev(cell, offset) += fun(i, qp) * value * dX;
                            }
                        });
                }
            }

            void describe(std::ostream &os) const override {
                const SizeType num_cells = fe_->num_cells();
                const int num_fields = fe_->num_fields();

                const int n_dofs = num_fields;

                std::cout << "num_cells: " << num_cells << ", num_fields: " << num_fields << "\n";

                for (SizeType c = 0; c < num_cells; ++c) {
                    os << c << ")\n";
                    for (SizeType i = 0; i < n_dofs; ++i) {
                        os << element_vectors_(c, i) << " ";

                        os << '\n';
                    }

                    os << '\n';
                }
            }

            inline const DynRankView &element_vectors() const { return element_vectors_; }

            // NVCC_PRIVATE :
            Op op_;
            std::shared_ptr<FE> fe_;
            DynRankView element_vectors_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FORCING_FUNCTION_HPP
