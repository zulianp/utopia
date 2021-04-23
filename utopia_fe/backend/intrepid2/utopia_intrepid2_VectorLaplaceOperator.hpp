#ifndef UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class Ceofficient>
    class VectorLaplaceOperator : public Configurable {
    public:
        void read(Input &in) override { in.get("coeff", coeff); }

        VectorLaplaceOperator(const Ceofficient &coeff = Ceofficient(1.0)) : coeff(coeff) {}
        Ceofficient coeff;
    };

    namespace intrepid2 {

        template <int Dim, typename Ceofficient, typename Scalar>
        class Assemble<VectorLaplaceOperator<Dim, Ceofficient>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::VectorLaplaceOperator<Dim, Ceofficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }
            int rank() const override { return 2; }
            inline std::string name() const override { return "VectorLaplaceOperator"; }

            bool assemble() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<VectorLaplaceOperator>::assemble");

                this->ensure_mat_accumulator();

                auto &fe = this->fe();

                // const int num_fields = fe.num_fields();
                // const int n_dofs = num_fields * fe.spatial_dimension();
                const int n_qp = fe.num_qp();

                auto data = this->data();

                // Only works if coeff is a scalar
                {
                    auto coeff = op_.coeff;
                    auto grad = fe.grad;
                    auto measure = fe.measure;

                    this->loop_cell_test_trial(
                        "Assemble<VectorLaplaceOperator>::init",
                        KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            StaticVector<Scalar, Dim> temp_i, temp_j;
                            StaticMatrix<Scalar, Dim, Dim> grad_i;
                            StaticMatrix<Scalar, Dim, Dim> grad_j;

                            // assert((Dim == 1 || &grad(cell, i, 0, 1) - &grad(cell, i, 0, 0) == 1UL) &&
                            //        "spatial dimension must be contiguos");


                            // grad: num_cells, n_fun, num_qp, spatial_dimension
                            // measure: num_cells, num_qp;

                            // needs loop over quadrature points and spatial_dim
                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto coeff_x_dX = coeff * measure(cell, qp);

                                for(int d = 0; d < Dim; ++d) {
                                    temp_i[d] = grad(cell, i, qp, d);
                                    temp_j[d] = grad(cell, j, qp, d);
                                }


                                for (int di = 0; di < Dim; ++di) {
                                    make_tensor_grad(di, &temp_i[0], grad_i);
                                    auto dof_i = i * Dim + di;

                                    for (int dj = 0; dj < Dim; ++dj) {
                                        auto dof_j = j * Dim + dj;

                                        make_tensor_grad(dj, &temp_j[0], grad_j);

                                        const Scalar val = inner(grad_i, grad_j) * coeff_x_dX;

                                        data(cell, dof_i, dof_j) += val;
                                    }
                                }
                            }
                        });
                }

                UTOPIA_TRACE_REGION_END("Assemble<VectorLaplaceOperator>::assemble");
                return true;
            }

            UTOPIA_INLINE_FUNCTION static void make_tensor_grad(const int dim,
                                                                Scalar *grad,
                                                                StaticMatrix<Scalar, Dim, Dim> &tensor_grad) {
                tensor_grad.set(0.0);

                for (int i = 0; i < Dim; ++i) {
                    tensor_grad(dim, i) = grad[i];
                }
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
