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

        Mass(const Fun &density) : density(density) {}

        UTOPIA_FUNCTION Mass() = default;
        UTOPIA_FUNCTION Mass(const Mass &) = default;

        Fun density;
        int n_components{1};
    };

    namespace intrepid2 {

        template <typename Fun>
        class Assemble<Mass<Fun>, typename Traits<Fun>::Scalar> : public Describable {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::Mass<Fun>;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : op_(std::move(op)), fe_(fe) {}

            void init() {
                const int n_components = op_.n_components;
                const int component = op_.component;

                const int num_fields = fe_->num_fields();
                const int n_dofs = num_fields * n_components;
                const int n_qp = fe_->num_qp();

                element_matrices_ = DynRankView("Mass", fe_->num_cells(), n_dofs, n_dofs);

                {
                    auto em = element_matrices_;
                    auto density = op_.density;
                    auto fun = fe_->fun;
                    auto measure = fe_->measure;

                    assert(n_components == 1 && "IMPLEMENT ME");

                    Kokkos::parallel_for(
                        "Assemble<Mass>::init",
                        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
                            {0, 0, 0}, {static_cast<int>(em.extent(0)), num_fields, num_fields}),
                        KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            // FIXME
                            auto offset_i = i * n_components;
                            auto offset_j = j * n_components;

                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);
                                em(cell, offset_i, offset_j) += fun(i, qp) * fun(j, qp) * density * dX;
                            }

                            // for(int sub_i = 0; sub_i < n_components; ++sub_i) {
                            //     for(int sub_j = 0; sub_j < n_components; ++sub_j) {
                            //         em(cell, offset_i + sub_i, offset_j + sub_j) = em(cell, offset_i, offset_j)
                            //     }
                            // }
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
                        for (SizeType j = 0; j < n_dofs; ++j) {
                            os << element_matrices_(c, i, j) << " ";
                        }

                        os << '\n';
                    }

                    os << '\n';
                }
            }

            inline const DynRankView &element_matrices() const { return element_matrices_; }

            // NVCC_PRIVATE :
            Op op_;
            std::shared_ptr<FE> fe_;
            DynRankView element_matrices_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID_2_MASS_HPP
