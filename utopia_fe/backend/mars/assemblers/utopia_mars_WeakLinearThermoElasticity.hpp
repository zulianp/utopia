#ifndef UTOPIA_MARS_WEAK_LINEAR_THERMO_ELASTICITY_HPP
#define UTOPIA_MARS_WEAK_LINEAR_THERMO_ELASTICITY_HPP

#include "utopia_Tracer.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"

#include "utopia_kokkos_LaplaceOp.hpp"
#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_kokkos_Strain.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"
#include "utopia_kokkos_SubView.hpp"

#include "utopia_kokkos_WeakLinearThermoElasticityCouplingOp.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class WeakLinearThermoElasticity : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;
            using FE = typename Super::FE;
            using Grad = typename FE::Gradient;
            using Measure = typename FE::Measure;
            using Function = typename FE::Function;
            static constexpr int Dim = DMesh::Dim;

            inline int n_vars() const { return Dim + 1; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    StressStrainParameters<Scalar, Scalar> ssp;
                    ssp.read(in);

                    lambda = ssp.first_lame_parameter.get();
                    mu = ssp.shear_modulus.get();

                    in.get("temperature", temperature);
                    in.get("displacement", displacement);

                    assert(temperature != displacement);
                }

                Params() {
                    const Scalar E = 50e3;
                    const Scalar nu = 0.2;

                    mu = E / 2 / (1 + nu);
                    lambda = E * nu / (1 + nu) / (1 - 2 * nu);
                    alpha = 1e-5;
                }

                Params(const Scalar &lambda, const Scalar &mu, const Scalar &alpha)
                    : lambda(lambda), mu(mu), alpha(alpha) {}

                Scalar lambda;
                Scalar mu;
                Scalar alpha;

                int displacement{0};
                int temperature{Dim};
            };

            using UOp = utopia::kokkos::kernels::LinearElasticityOp<Dim, Scalar, Scalar, Scalar, Grad, Measure>;

            using COp = utopia::kokkos::kernels::LaplaceOp<Scalar, Scalar, Grad, Measure>;

            using UCOp =
                utopia::kokkos::kernels::WeakLinearThermoElasticityCouplingOp<Dim, Scalar, Grad, Function, Measure>;

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool is_operator() const override { return true; }

            bool apply(const Vector &x, Vector &vec) override {
                this->ensure_fe();

                auto &fe = this->fe();

                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int displacement = params_.displacement;
                const int temperature = params_.temperature;
                const int n_var = this->n_vars();

                UOp u_op(params_.lambda, params_.mu, grad, measure);
                COp c_op(1., grad, measure);
                UCOp uc_op(params_.lambda, params_.mu, params_.alpha, grad, fun, measure);

                return Super::coupled_vector_scalar_op_apply(displacement, temperature, u_op, c_op, uc_op, x, vec);
            }

            bool assemble(Matrix &hessian) override {
                this->ensure_fe();

                auto &fe = this->fe();

                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int displacement = params_.displacement;
                const int temperature = params_.temperature;
                const int n_var = this->n_vars();

                UOp u_op(params_.lambda, params_.mu, grad, measure);
                COp c_op(1., grad, measure);
                UCOp uc_op(params_.lambda, params_.mu, params_.alpha, grad, fun, measure);

                return Super::coupled_vector_scalar_op_assemble(displacement, temperature, u_op, c_op, uc_op, hessian);
            }

            void read(Input &in) override { params_.read(in); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) override {
                Super::set_environment(env);
                // DO extra things ...
            }

            inline bool is_linear() const override { return true; }

        private:
            Params params_;
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_WEAK_LINEAR_THERMO_ELASTICITY_HPP