#ifndef UTOPIA_MARS_LINEAR_ELASTICITY_HPP
#define UTOPIA_MARS_LINEAR_ELASTICITY_HPP

#include "mars_quad4.hpp"
#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_mars_ConcreteFEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class LinearElasticity : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;
            static constexpr int Dim = DMesh::Dim;

            // Correct fe type is managed in ConcreteFEAssembler
            using FE = typename Super::FE;
            using Op = utopia::kokkos::kernels::
                LinearElasticityOp<Dim, Scalar, Scalar, Scalar, typename FE::Gradient, typename FE::Measure>;

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool is_operator() const override { return true; }

            bool apply(const Vector &x, Vector &vec) override {
                this->ensure_fe();

                auto &fe = this->fe();

                // This is for linear materials since gradient := A x - b
                // non linear materials will have to implement a separate gradient Op
                Op op(lambda_, mu_, fe.grad(), fe.measure());
                return Super::block_op_apply(op, x, vec);
            }

            bool assemble(Matrix &hessian) override {
                this->ensure_fe();

                auto &fe = this->fe();
                Op op(lambda_, mu_, fe.grad(), fe.measure());
                return Super::block_op_assemble_matrix(op, hessian);
            }

            void read(Input &in) override {
                in.get("lambda", lambda_);
                in.get("mu", mu_);
            }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) override {
                Super::set_environment(env);
                // DO extra things ...
            }

            inline bool is_linear() const override { return true; }

        private:
            Scalar lambda_{1};
            Scalar mu_{1};
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LINEAR_ELASTICITY_HPP
