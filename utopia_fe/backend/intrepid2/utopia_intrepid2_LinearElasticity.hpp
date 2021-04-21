#ifndef UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
#define UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class LinearElasticity : public Configurable {
    public:
        void read(Input &in) override {
            in.get("lambda", lambda);
            in.get("mu", mu);
        }

        LinearElasticity(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                         const ShearModulus &mu = FirstLameParameter(1.0))
            : lambda(lambda), mu(mu) {}

        FirstLameParameter lambda;
        ShearModulus mu;
    };

    namespace intrepid2 {

        template <int Dim, class FirstLameParameter, class ShearModulus, typename Scalar>
        class Assemble<LinearElasticity<Dim, FirstLameParameter, ShearModulus>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::LinearElasticity<Dim, FirstLameParameter, ShearModulus>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            int rank() const override { return 2; }
            inline std::string name() const override { return "LinearElasticity"; }
            bool assemble() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LinearElasticity>::assemble");

                this->ensure_mat_accumulator();

                auto &fe = this->fe();
                auto data = this->data();

                // const int num_fields = fe_.num_fields();
                // const int n_dofs = num_fields * fe_.spatial_dimension();
                const int n_qp = fe.num_qp();

                // Only works if coeff is a scalar
                {
                    auto mu = op_.mu;
                    auto lambda = op_.lambda;
                    auto grad = fe.grad;
                    auto mux2 = mu * 2;
                    auto measure = fe.measure;

                    this->loop_cell_test_trial(
                        "Assemble<LinearElasticity>::init", KOKKOS_LAMBDA(const int &cell, const int &i, const int &j) {
                            StaticMatrix<Scalar, Dim, Dim> strain_i;
                            StaticMatrix<Scalar, Dim, Dim> strain_j;

                            assert((Dim == 1 || &grad(cell, i, 0, 1) - &grad(cell, i, 0, 0) == 1UL) &&
                                   "spatial dimension must be contiguos");

                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);

                                for (int di = 0; di < Dim; ++di) {
                                    make_strain(di, &grad(cell, i, qp, 0), strain_i);
                                    const Scalar trace_i = trace(strain_i);
                                    auto dof_i = i * Dim + di;

                                    for (int dj = 0; dj < Dim; ++dj) {
                                        auto dof_j = j * Dim + dj;

                                        make_strain(dj, &grad(cell, j, qp, 0), strain_j);

                                        const Scalar val =
                                            (mux2 * inner(strain_i, strain_j) + lambda * trace_i * trace(strain_j)) *
                                            dX;

                                        data(cell, dof_i, dof_j) += val;
                                    }
                                }
                            }
                        });
                }

                UTOPIA_TRACE_REGION_END("Assemble<LinearElasticity>::assemble");
                return true;
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

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
