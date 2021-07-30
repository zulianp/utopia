#ifndef UTOPIA_MARS_LAPLACE_OPERATOR_HPP
#define UTOPIA_MARS_LAPLACE_OPERATOR_HPP

#include "mars_quad4.hpp"
#include "utopia_kokkos_LaplaceOp.hpp"
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

            // Correct fe type is managed in ConcreteFEAssembler
            using FE = typename Super::FE;

            // Same op template used in intrepid2-based backend
            using Op = utopia::kokkos::kernels::LaplaceOp<Scalar, Scalar, typename FE::Gradient, typename FE::Measure>;

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                if (!assemble(hessian)) {
                    return false;
                }

                assemble(x, gradient);
                return true;
            }

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool assemble(const Vector &x, Vector &vec) override {
                this->ensure_fe();

                auto &fe = this->fe();

                // This is for linear materials since gradient := A x - b
                // non linear materials will have to implement a separate gradient Op
                Op op(coeff_, fe.grad(), fe.measure());
                return Super::scalar_op_apply(op, x, vec);
            }

            bool assemble(Matrix &hessian) override {
                this->ensure_fe();

                auto &fe = this->fe();
                Op op(coeff_, fe.grad(), fe.measure());
                return Super::scalar_op_assemble_matrix(op, hessian);
            }

            void read(Input &in) override { in.get("coeff", coeff_); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) override {
                Super::set_environment(env);
                // DO extra things ...
            }

            inline bool is_linear() const override { return true; }

        private:
            Scalar coeff_{1.0};
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP
