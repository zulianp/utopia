#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_kokkos_LaplaceOp.hpp"

namespace utopia {

    namespace kokkos {

        template <class FunctionSpace, class FE_, class DiffusionCoefficient = typename FE_::Scalar>
        class LaplaceOperator : public FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>>;
            using Gradient = typename FE::Gradient;
            using Measure = typename FE::Measure;

            using Params = utopia::kokkos::LaplaceOp<FE, DiffusionCoefficient>;
            using Op = typename Params::UniformKernel;

            LaplaceOperator(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "LaplaceOperator"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline Op make_op() { return op_.uniform_kernel(this->fe()); }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LaplaceOperator>::assemble");
                this->ensure_matrix_accumulator();

                this->loop_cell_test_trial("Assemble<LaplaceOperator>::assemble",
                                           op_and_store_cell_ij(this->matrix_data(), op_.uniform_kernel(this->fe())));

                if (!op_.subdomain_value.empty) {
                    auto data = this->matrix_data();
                    this->scale_matrix(op_.subdomain_value, data);
                }

                assert(check_op());

                UTOPIA_TRACE_REGION_END("Assemble<LaplaceOperator>::assemble");
                return true;
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LaplaceOperator>::apply");

                this->apply_operator("Assemble<LaplaceOperator>::apply", x, y, make_op());

                if (!op_.subdomain_value.empty) {
                    auto data = this->vector_data();
                    this->scale_vector(op_.subdomain_value, data);
                }

                UTOPIA_TRACE_REGION_END("Assemble<LaplaceOperator>::apply");
                return true;
            }

            bool check_op() {
                auto data = this->matrix_data();
                const int n_shape_functions = this->fe().n_shape_functions();

                this->loop_cell(
                    "LaplaceOperator::check_op()", UTOPIA_LAMBDA(int cell) {
                        for (int i = 0; i < n_shape_functions; ++i) {
                            for (int j = 0; j < n_shape_functions; ++j) {
                                auto val = data(cell, i, j);
                                assert(i != j || val > 0.0);
                            }
                        }
                    });

                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif
