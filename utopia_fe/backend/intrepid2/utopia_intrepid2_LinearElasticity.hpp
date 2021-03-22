#ifndef UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
#define UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP

#include "utopia_intrepid2_LaplaceOperator.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class LinearElasticity {
    public:
        FirstLameParameter lambda{1.0};
        ShearModulus mu{1.0};
    };

    namespace intrepid2 {

        template <int Dim, class FirstLameParameter, class ShearModulus, typename Scalar>
        class Assemble<LinearElasticity<Dim, FirstLameParameter, ShearModulus>, Scalar> : public Describable {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using Op = utopia::LinearElasticity<Dim, FirstLameParameter, ShearModulus>;
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
                    auto mu = op_.mu;
                    auto lambda = op_.lambda;
                    auto grad = fe_->grad;
                    auto mux2 = mu * 2;
                    auto measure = fe_->measure;

                    Kokkos::parallel_for(
                        "scale_with_coeff",
                        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
                            {0, 0, 0}, {static_cast<int>(em.extent(0)), num_fields, num_fields}),
                        KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            StaticMatrix<Scalar, Dim, Dim> strain_i;
                            StaticMatrix<Scalar, Dim, Dim> strain_j;

                            assert(&grad(cell, i, j, 1) - &grad(cell, i, j, 0) == 1UL &&
                                   "spatial dimension must be contiguos");

                            // grad: num_cells, n_fun, num_qp, spatial_dimension
                            // measure: num_cells, num_qp;

                            // needs loop over quadrature points and spatial_dim
                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);

                                for (int di = 0; di < Dim; ++di) {
                                    make_strain(di, &grad(cell, i, qp), strain_i);
                                    const Scalar trace_i = trace(strain_i);
                                    auto dof_i = i * Dim + di;

                                    em(cell, dof_i, dof_i) +=
                                        (mux2 * inner(strain_i, strain_i) + lambda * trace_i * trace_i) * dX;

                                    for (int dj = di + 1; dj < Dim; ++dj) {
                                        auto dof_j = j * Dim + dj;

                                        make_strain(dj, &grad(cell, j, qp), strain_j);

                                        const Scalar val =
                                            (mux2 * inner(strain_i, strain_j) + lambda * trace_i * trace(strain_j)) *
                                            dX;

                                        em(cell, dof_i, dof_j) += val;
                                        em(cell, dof_j, dof_i) += val;
                                    }
                                }
                            }
                        });
                }
            }

            UTOPIA_INLINE_FUNCTION static void make_strain(const int dim,
                                                           Scalar *grad,
                                                           StaticMatrix<Scalar, Dim, Dim> &strain) {
                strain.set(0.0);

                for (int i = 0; i < Dim; ++i) {
                    strain(dim, i) = grad[i];
                }

                strain.symmetrize();
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

#endif  // UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
