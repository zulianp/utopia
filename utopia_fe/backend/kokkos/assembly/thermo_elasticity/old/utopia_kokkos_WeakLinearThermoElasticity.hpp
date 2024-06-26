#ifndef UTOPIA_KOKKOS_WEAK_LINEAR_THERMO_ELASTICITY_HPP
#define UTOPIA_KOKKOS_WEAK_LINEAR_THERMO_ELASTICITY_HPP

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

    namespace kokkos {

        template <class FunctionSpace,
                  class FE_,
                  int Dim,
                  class FirstLameParameter = typename FE_::Scalar,
                  class ShearModulus = FirstLameParameter,
                  class Alpha = FirstLameParameter>
        class WeakLinearThermoElasticity
            : public utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using Grad = typename FE::Gradient;
            using Measure = typename FE::Measure;
            using Function = typename FE::Function;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>>;
            using VectorView = typename Super::VectorView;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
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

                Params(const FirstLameParameter &lambda, const ShearModulus &mu, const Alpha &alpha)
                    : lambda(lambda), mu(mu), alpha(alpha) {}

                FirstLameParameter lambda;
                ShearModulus mu;
                Alpha alpha;

                int displacement{0};
                int temperature{Dim};
            };

            using UOp = utopia::kokkos::kernels::
                LinearElasticityOp<Dim, Scalar, FirstLameParameter, ShearModulus, Grad, Measure>;

            using COp = utopia::kokkos::kernels::LaplaceOp<Scalar, Scalar, Grad, Measure>;

            using UCOp =
                utopia::kokkos::kernels::WeakLinearThermoElasticityCouplingOp<Dim, Scalar, Grad, Function, Measure>;

            // class UCOp {
            // public:
            //     using StrainKernel = utopia::kokkos::kernels::LinearizedStrain<Dim, Scalar>;

            //     UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
            //                                              const int i,
            //                                              const int j,
            //                                              const int sub_i) const {
            //         Scalar ret = 0.0;
            //         for (int qp = 0; qp < n_quad_points; ++qp) {
            //             ret += (*this)(cell, i, j, qp, sub_i);
            //         }
            //         return ret;
            //     }

            //     UTOPIA_INLINE_FUNCTION Scalar
            //     operator()(const int cell, const int i, const int j, const int qp, const int sub_i) const {
            //         auto v = -alpha * (3 * lambda + 2 * mu) * fun(j, qp);
            //         return StrainKernel::trace(grad, cell, i, qp, sub_i) * v * measure(cell, qp);
            //     }

            //     UTOPIA_INLINE_FUNCTION UCOp(const Scalar &lambda,
            //                                 const Scalar &mu,
            //                                 const Scalar &alpha,
            //                                 const Grad &grad,
            //                                 const Function &fun,
            //                                 const Measure &measure)
            //         : lambda(lambda),
            //           mu(mu),
            //           alpha(alpha),
            //           grad(grad),
            //           fun(fun),
            //           measure(measure),
            //           n_quad_points(measure.extent(1)) {}

            //     Scalar lambda, mu, alpha;
            //     Grad grad;
            //     Function fun;
            //     Measure measure;
            //     int n_quad_points;
            // };

            WeakLinearThermoElasticity(const std::shared_ptr<FE> &fe, Params op = Params())
                : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim + 1; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline std::string name() const override { return "WeakLinearThermoElasticity"; }

            bool apply(const VectorView &x, VectorView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("WeakLinearThermoElasticity::apply");

                auto &fe = this->fe();
                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int displacement = op_.displacement;
                const int temperature = op_.temperature;
                const int n_var = this->n_vars();

                UOp u_op(op_.lambda, op_.mu, grad, measure);
                COp c_op(1., grad, measure);
                UCOp uc_op(op_.lambda, op_.mu, op_.alpha, grad, fun, measure);

                this->loop_cell_test(
                    "WeakLinearThermoElasticity::apply", UTOPIA_LAMBDA(const int &cell, const int &i) {
                        for (int j = 0; j < n_shape_functions; ++j) {
                            for (int sub_i = 0; sub_i < Dim; ++sub_i) {
                                Scalar val = 0.0;

                                // Integrate displacement block
                                for (int sub_j = 0; sub_j < Dim; ++sub_j) {
                                    const int idx_j = j * n_var + displacement + sub_j;
                                    val += u_op(cell, i, j, sub_i, sub_j) * x(cell, idx_j);
                                }

                                // Integrate coupling block
                                val += uc_op(cell, i, j, sub_i) * x(cell, j * n_var + temperature);

                                // Store result in displacement block
                                y(cell, i * n_var + displacement + sub_i) += val;
                            }

                            // Integrate and store in temperature block
                            y(cell, i * n_var + temperature) += c_op(cell, i, j) * x(cell, j * n_var + temperature);
                        }
                    });

                UTOPIA_TRACE_REGION_END("WeakLinearThermoElasticity::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("WeakLinearThermoElasticity::assemble_matrix");
                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int displacement = op_.displacement;
                const int temperature = op_.temperature;
                const int n_var = this->n_vars();

                UOp u_op(op_.lambda, op_.mu, grad, measure);
                COp c_op(1., grad, measure);
                UCOp uc_op(op_.lambda, op_.mu, op_.alpha, grad, fun, measure);

                auto data = this->matrix_data();

                this->loop_cell_test_trial(
                    "WeakLinearThermoElasticity::assemble_matrix",
                    UTOPIA_LAMBDA(const int &cell, const int &i, const int &j) {
                        const int offset_i = i * n_var;
                        const int offset_j = j * n_var;

                        for (int sub_i = 0; sub_i < Dim; ++sub_i) {
                            int u_dof_i = offset_i + sub_i + displacement;

                            // Integrate displacement block
                            for (int sub_j = 0; sub_j < Dim; ++sub_j) {
                                int dof_j = offset_j + sub_j + displacement;
                                data(cell, u_dof_i, dof_j) = u_op(cell, i, j, sub_i, sub_j);
                            }

                            // Integrate coupling block
                            int dof_j = offset_j + temperature;
                            data(cell, u_dof_i, dof_j) = uc_op(cell, i, j, sub_i);

                            // Store result in displacement block
                            data(cell, offset_i + temperature, dof_j) += c_op(cell, i, j);
                        }
                    });

                UTOPIA_TRACE_REGION_END("WeakLinearThermoElasticity::assemble_matrix");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_WEAK_LINEAR_THERMO_ELASTICITY_HPP