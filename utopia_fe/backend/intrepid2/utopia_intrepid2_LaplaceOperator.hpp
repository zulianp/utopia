#ifndef UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_LAPLACE_OPERATOR_HPP

#include "utopia_fe_base.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"

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

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "LaplaceOperator"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const DiffusionCoefficient &coeff,
                                          const DynRankView &grad,
                                          const DynRankView &measure)
                    : coeff(coeff), grad(grad), measure(measure), n_qp(measure.extent(1)), dim(grad.extent(3)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar dot_g = 0.0;
                        for (int d = 0; d < dim; ++d) {
                            dot_g += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += dot_g * measure(cell, qp);
                    }

                    return ret * coeff;
                }

                const DiffusionCoefficient coeff;
                const DynRankView grad;
                const DynRankView measure;
                const int n_qp;
                const int dim;
            };

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.coeff, fe.grad, fe.measure);
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

            // bool assemble_matrix_intrepid_tutorial() {
            //     UTOPIA_TRACE_REGION_BEGIN("Assemble<LaplaceOperator>::assemble_matrix_intrepid_tutorial");
            //     this->ensure_matrix_accumulator();

            //     auto &fe = this->fe();
            //     auto data = this->matrix_data();

            //     DynRankView grad_x_measure(
            //         "grad_x_measure", fe.num_cells(), fe.num_fields(), fe.num_qp(), fe.spatial_dimension());

            //     Kokkos::deep_copy(grad_x_measure, fe.grad);
            //     FunctionSpaceTools::template multiplyMeasure<Scalar>(grad_x_measure, fe.measure, fe.grad);
            //     FunctionSpaceTools::template integrate<Scalar>(data, fe.grad, grad_x_measure);

            //     // Only works if coeff is a scalar
            //     if (op_.coeff != 1.0) {
            //         auto c = op_.coeff;
            //         this->loop_cell_test_trial(
            //             "scale_with_coeff",
            //             KOKKOS_LAMBDA(const int &i0, const int &i1, const int &i2) { data(i0, i1, i2) *= c; });
            //     }

            //     UTOPIA_TRACE_REGION_END("Assemble<LaplaceOperator>::assemble_matrix_intrepid_tutorial");
            //     return true;
            // }

            // NVCC_PRIVATE :
            UserOp op_;
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif
