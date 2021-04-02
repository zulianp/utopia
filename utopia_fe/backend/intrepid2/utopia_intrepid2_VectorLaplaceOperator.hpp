#ifndef UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
#define UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP

#include "utopia_intrepid2_LaplaceOperator.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class Ceofficient>
    class VectorLaplaceOperator : public Configurable {
    public:
        void read(Input &in) override { in.get("coeff", coeff); }

        VectorLaplaceOperator(const Ceofficient &coeff) : coeff(coeff) {}
        Ceofficient coeff;
    };

    namespace intrepid2 {

        template <int Dim, typename Ceofficient, typename Scalar>
        class Assemble<VectorLaplaceOperator<Dim, Ceofficient>, Scalar> : public Describable {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::VectorLaplaceOperator<Dim, Ceofficient>;
            using ExecutionSpace = typename FE::ExecutionSpace;

            Assemble(Op op, const std::shared_ptr<FE> &fe) : op_(std::move(op)), fe_(fe) {}

            void init() {
                assert(Dim == fe_->spatial_dimension());
                const int num_fields = fe_->num_fields();
                const int n_dofs = num_fields * fe_->spatial_dimension();
                const int n_qp = fe_->num_qp();

                element_matrices_ = DynRankView("stiffness_matrix", fe_->num_cells(), n_dofs, n_dofs);

                // Only works if coeff is a scalar
                {
                    auto em = element_matrices_;
                    auto coeff = op_.coeff;
                    auto grad = fe_->grad;
                    auto measure = fe_->measure;

                    Kokkos::parallel_for(
                        "Assemble<VectorLaplaceOperator>::init",
                        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
                            {0, 0, 0}, {static_cast<int>(em.extent(0)), num_fields, num_fields}),
                        KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            StaticMatrix<Scalar, Dim, Dim> grad_i;
                            StaticMatrix<Scalar, Dim, Dim> grad_j;

                            assert((Dim == 1 || &grad(cell, i, 0, 1) - &grad(cell, i, 0, 0) == 1UL) &&
                                   "spatial dimension must be contiguos");

                            // grad: num_cells, n_fun, num_qp, spatial_dimension
                            // measure: num_cells, num_qp;

                            // needs loop over quadrature points and spatial_dim
                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto coeff_x_dX = coeff * measure(cell, qp);

                                for (int di = 0; di < Dim; ++di) {
                                    make_tensor_grad(di, &grad(cell, i, qp, 0), grad_i);
                                    auto dof_i = i * Dim + di;

                                    for (int dj = 0; dj < Dim; ++dj) {
                                        auto dof_j = j * Dim + dj;

                                        make_tensor_grad(dj, &grad(cell, j, qp, 0), grad_j);

                                        const Scalar val = inner(grad_i, grad_j) * coeff_x_dX;

                                        em(cell, dof_i, dof_j) += val;
                                    }
                                }
                            }
                        });
                }
            }

            UTOPIA_INLINE_FUNCTION static void make_tensor_grad(const int dim,
                                                                Scalar *grad,
                                                                StaticMatrix<Scalar, Dim, Dim> &tensor_grad) {
                tensor_grad.set(0.0);

                for (int i = 0; i < Dim; ++i) {
                    tensor_grad(dim, i) = grad[i];
                }
            }

            void describe(std::ostream &os) const override {
                const SizeType num_cells = fe_->num_cells();
                const int num_fields = fe_->num_fields();

                const int n_dofs = num_fields * fe_->spatial_dimension();

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

#endif  // UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
