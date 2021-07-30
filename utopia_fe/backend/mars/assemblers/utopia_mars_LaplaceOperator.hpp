#ifndef UTOPIA_MARS_LAPLACE_OPERATOR_HPP
#define UTOPIA_MARS_LAPLACE_OPERATOR_HPP

#include "mars_quad4.hpp"
#include "utopia_kokkos_LaplaceOp.hpp"
#include "utopia_mars_ConcreteFEAssembler.hpp"
#include "utopia_mars_FEBuilder.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class LaplaceOperator : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;

            using FEHandler = typename Super::FEHandler;
            using DofHandler = typename Super::DofHandler;
            using FEDofMap = typename Super::FEDofMap;
            using SPattern = typename Super::SPattern;

            static const ::mars::Integer Type = SPattern::DofHandler::ElemType;
            using UniformFE = utopia::kokkos::UniformFE<Scalar>;

            // Same op template used in intrepid2-based backend
            using Op = utopia::kokkos::kernels::
                LaplaceOp<Scalar, Scalar, typename UniformFE::Gradient, typename UniformFE::Measure>;

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                if (!assemble(hessian)) {
                    return false;
                }

                gradient = hessian * x;
                return true;
            }

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool assemble(const Vector &x, Vector &vec) override {
                assert(false && "IMPLEMENT ME");
                return false;
            }

            bool assemble(Matrix &hessian) override {
                init();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                Op op(coeff_, fe_->grad(), fe_->measure());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, i);
                        if (!dof_handler.is_owned(local_dof_i)) continue;

                        for (int j = 0; j < n_fun; j++) {
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, j);
                            Scalar val = op(elem_index, i, j);

                            assert(val == val);
                            sp.atomic_add_value(local_dof_i, local_dof_j, val);
                        }
                    }
                });

                // // FIXME Bad it should not do this
                auto mat_impl =
                    Teuchos::rcp(new Matrix::CrsMatrixType(this->handler()->get_sparsity_pattern().get_matrix(),
                                                           hessian.raw_type()->getRowMap(),
                                                           hessian.raw_type()->getColMap()));
                hessian.wrap(mat_impl, false);
                return true;
            }

            void read(Input &in) override { in.get("coeff", coeff_); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &) override {
                // assert(false);
            }

            inline bool is_linear() const override { return true; }

        private:
            Scalar coeff_{1.0};
            std::unique_ptr<UniformFE> fe_;

            void init() {
                if (!fe_) {
                    fe_ = utopia::make_unique<UniformFE>();
                    FEBuilder<FEHandler, utopia::kokkos::UniformFE<Scalar>> builder;
                    builder.build(*this->handler(), *fe_);
                }
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP
