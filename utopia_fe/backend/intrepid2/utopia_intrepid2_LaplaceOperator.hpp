#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"

#include "utopia_kokkos_LaplaceOp.hpp"

namespace utopia {

    template <class DiffusionCoefficient>
    class LaplaceOperator : public Configurable {
    public:
        void read(Input &in) override {
            in.get("coeff", coeff);
            in.get("permeability_function", [this](Input &node) {
                node.get("function", [this](Input &function_node) { subdomain_function.read(function_node); });
            });

            in.get("diffusion_function", [this](Input &node) {
                node.get("function", [this](Input &function_node) { subdomain_function.read(function_node); });
            });

            bool verbose = false;
            in.get("verbose", verbose);

            if (verbose) {
                utopia::out() << "------------------------------\n";
                utopia::out() << "LaplaceOperator:\n";
                utopia::out() << "coeff:\t" << coeff << '\n';

                if (!subdomain_function.empty) {
                    utopia::out() << "function: SubdomainValue\n";
                }

                utopia::out() << "------------------------------\n";
            }
        }

        LaplaceOperator(const DiffusionCoefficient &coeff = DiffusionCoefficient(1.0))
            : coeff(coeff), subdomain_function(1.0) {}
        DiffusionCoefficient coeff;

        // FIXME
        intrepid2::SubdomainValue<DiffusionCoefficient> subdomain_function;
    };

    namespace intrepid2 {

        template <class DiffusionCoefficient, typename Scalar>
        class Assemble<LaplaceOperator<DiffusionCoefficient>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::LaplaceOperator<DiffusionCoefficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;
            using Gradient = typename FE::Gradient;
            using Measure = typename FE::Measure;

            using Op = utopia::kokkos::kernels::LaplaceOp<Scalar, DiffusionCoefficient, Gradient, Measure>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "LaplaceOperator"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.coeff, fe.grad(), fe.measure());
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LaplaceOperator>::assemble");
                this->ensure_matrix_accumulator();

                this->loop_cell_test_trial("Assemble<LaplaceOperator>::assemble",
                                           op_and_store_cell_ij(this->matrix_data(), make_op()));

                if (!op_.subdomain_function.empty) {
                    auto data = this->matrix_data();
                    this->scale_matrix(op_.subdomain_function, data);
                }

                assert(check_op());

                UTOPIA_TRACE_REGION_END("Assemble<LaplaceOperator>::assemble");
                return true;

                // return assemble_matrix_intrepid_tutorial();
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LaplaceOperator>::apply");

                this->apply_operator("Assemble<LaplaceOperator>::apply", x, y, make_op());

                if (!op_.subdomain_function.empty) {
                    auto data = this->vector_data();
                    this->scale_vector(op_.subdomain_function, data);
                }

                UTOPIA_TRACE_REGION_END("Assemble<LaplaceOperator>::apply");
                return true;
            }

            bool check_op() {
                auto data = this->matrix_data();
                const int n_shape_functions = this->fe().n_shape_functions();

                // this->matrix_accumulator()->describe(utopia::out().stream());

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
            UserOp op_;
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif
