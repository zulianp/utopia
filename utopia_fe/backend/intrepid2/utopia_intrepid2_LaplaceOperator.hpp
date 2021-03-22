#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include "utopia_intrepid2_FE.hpp"

#include <Kokkos_DynRankView.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_Cubature.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>

namespace utopia {

    template <class DiffusionCoefficient>
    class LaplaceOperator {
    public:
        DiffusionCoefficient coeff;
    };

    namespace intrepid2 {

        template <class Operator, typename Scalar = UScalar>
        class Assemble {};

        template <class DiffusionCoefficient, typename Scalar>
        class Assemble<LaplaceOperator<DiffusionCoefficient>, Scalar> : public Describable {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::LaplaceOperator<DiffusionCoefficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : op_(std::move(op)), fe_(fe) {}

            void init() {
                element_matrices_ =
                    DynRankView("stiffness_matrix", fe_->num_cells(), fe_->num_fields(), fe_->num_fields());

                DynRankView grad_x_measure(
                    "grad_x_measure", fe_->num_cells(), fe_->num_fields(), fe_->num_qp(), fe_->spatial_dimension());

                Kokkos::deep_copy(grad_x_measure, fe_->grad);
                FunctionSpaceTools::template multiplyMeasure<Scalar>(grad_x_measure, fe_->measure, fe_->grad);
                FunctionSpaceTools::template integrate<Scalar>(element_matrices_, fe_->grad, grad_x_measure);

                // Only works if coeff is a scalar
                {
                    auto em = element_matrices_;
                    auto c = op_.coeff;

                    Kokkos::parallel_for(
                        "scale_with_coeff",
                        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
                            {0, 0, 0}, {em.extent(0), em.extent(1), em.extent(2)}),
                        KOKKOS_LAMBDA(const int &i0, const int &i1, const int &i2) { em(i0, i1, i2) *= c; });
                }
            }

            void describe(std::ostream &os) const override {
                const SizeType num_cells = fe_->num_cells();
                const SizeType num_fields = fe_->num_fields();

                std::cout << "num_cells: " << num_cells << ", num_fields: " << num_fields << "\n";

                for (SizeType c = 0; c < num_cells; ++c) {
                    os << c << ")\n";
                    for (SizeType i = 0; i < num_fields; ++i) {
                        for (SizeType j = 0; j < num_fields; ++j) {
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

#endif