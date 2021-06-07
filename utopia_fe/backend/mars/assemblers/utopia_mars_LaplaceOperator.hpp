#ifndef UTOPIA_MARS_LAPLACE_OPERATOR_HPP
#define UTOPIA_MARS_LAPLACE_OPERATOR_HPP

#include "mars_quad4.hpp"
#include "utopia_mars_ConcreteFEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class LaplaceOperator : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;
            using FE = typename Super::FE;

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
                fe_->intrepid2_laplace_operator();

                // FIXME Bad it should not do this
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
            std::unique_ptr<FE> fe_;

            void init() {
                if (!fe_) {
                    auto handler = this->handler();
                    fe_ = utopia::make_unique<FE>(handler->get_sparsity_pattern(), handler->get_fe_dof_map());
                    fe_->init();
                }
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP
