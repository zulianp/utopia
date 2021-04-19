#ifndef UTOPIA_INTREPID2_GRADIENT_HPP
#define UTOPIA_INTREPID2_GRADIENT_HPP

#include "utopia_Traits.hpp"
#include "utopia_Views.hpp"

#include "utopia_intrepid2_Commons.hpp"
#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

namespace utopia {
    namespace intrepid2 {

        template <typename Scalar>
        class Gradient {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Gradient(const std::shared_ptr<FE> &fe) : fe_(fe) {}

            class Op {
            public:
                Op(const DynRankView &grad, const DynRankView &coeff, DynRankView &field)
                    : grad(grad), coeff(coeff), field(field), num_fields(grad.extent(1)) {}

                void operator()(const int cell, const int qp, const int d) const {
                    for (int i = 0; i < num_fields; ++i) {
                        field(cell, qp, d) += coeff(cell, i) * grad(cell, i, qp, d);
                    }
                }

                const DynRankView grad, coeff;
                DynRankView field;
                const int num_fields;
            };

            void init(const DynRankView &coeff) {
                ensure_field();
                Kokkos::parallel_for(name(), grad_range(), Op(fe_->grad, coeff, field_));
            }

            inline std::string name() const { return "Gradient"; }

            void ensure_field() {
                if (field_.extent(0) < fe_->num_cells() || field_.extent(1) < fe_->num_qp()) {
                    field_ = DynRankView("gradient", fe_->num_cells(), fe_->num_qp());
                } else {
                    fill(field_, 0.0);
                }
            }

            inline Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace> grad_range() {
                int num_cells = fe_->num_cells();
                int num_qp = fe_->num_qp();
                int spatial_dimension = fe_->spatial_dimension();

                return Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0},
                                                                              {num_cells, num_qp, spatial_dimension});
            }

            inline DynRankView &field() { return field_; }
            inline std::shared_ptr<FE> fe() { return fe_; }

        private:
            std::shared_ptr<FE> fe_;
            DynRankView field_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_GRADIENT_HPP
