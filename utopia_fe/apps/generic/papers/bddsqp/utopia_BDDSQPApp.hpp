#ifndef UTOPIA_BDDSQP_APP_HPP
#define UTOPIA_BDDSQP_APP_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Field.hpp"
#include "utopia_MPITimeStatistics.hpp"

#include "utopia_ui.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include "utopia_fe_Core.hpp"

#include "utopia_Material.hpp"
#include "utopia_kokkos_MaterialFactory_impl.hpp"

#include <memory>

namespace utopia {

    template <class FunctionSpace, class FE>
    class BDDSQPApp {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using MaterialFactory_t = utopia::kokkos::MaterialFactory<FunctionSpace, FE>;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        using NewtonBase_t = utopia::NewtonBase<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        static void print_norms(const Matrix_t &hessian, const Vector_t &grad, const Vector_t &x) {
            Scalar_t norm_x = norm2(x);
            Scalar_t norm_g = norm2(grad);
            Scalar_t sum_g = sum(grad);
            Scalar_t norm_h = norm2(hessian);
            Scalar_t sum_h = sum(hessian);

            if (!x.comm().rank()) {
                utopia::out() << "-------------------------------\n";
                utopia::out() << "norm_h  : " << norm_h << "\n";
                utopia::out() << "sum_h   : " << sum_h << "\n";
                utopia::out() << "norm_g  : " << norm_g << "\n";
                utopia::out() << "norm_x  : " << norm_x << "\n";
                utopia::out() << "sum_g   : " << sum_g << "\n";
                utopia::out() << "-------------------------------\n";
            }
        }

        void run(Input &in) {
            FunctionSpace space;
            Field<FunctionSpace> input_state;
            bool use_state_as_initial_condition = false;

            in.get("use_state_as_initial_condition", use_state_as_initial_condition);

            in.require("space", [&](Input &in) {
                // bool read_state = false;
                // in.get("read_state", read_state);

                // if (read_state) {
                //     space.read_with_state(in, input_state);

                //     const Scalar_t norm_input_state = norm2(input_state.data());
                //     utopia::out() << "norm_input_state: " << norm_input_state << std::endl;

                // } else {
                space.read(in);
                // }
            });

            if (space.empty()) {
                Utopia::Abort("Space is empty!");
                return;
            }

            std::unique_ptr<AbstractMaterial<FunctionSpace>> material;

            in.require("problem", [&](Input &problem_node) {
                problem_node.require("assembly", [&](Input &assembly_node) {
                    assembly_node.require("material", [&](Input &material_node) {
                        std::string material_type;
                        material_node.require("type", material_type);
                        material = MaterialFactory_t::make(space.mesh().spatial_dimension(), material_type);
                        material->initialize(make_ref(space));

                        if (!material) {
                            Utopia::Abort();
                        }

                        material->read(material_node);
                    });
                });
            });

            Path output = "out.e";
            in.get("output", output);

            bool debug = false;
            in.get("debug", debug);

            ////////////////////////////////////////////////////////////////////////////////////////////////

            auto cg = std::make_shared<ConjugateGradient<Matrix_t, Vector_t, HOMEMADE>>();
            cg->apply_gradient_descent_step(true);
            std::shared_ptr<LinearSolver_t> linear_solver = cg;

            // std::shared_ptr<LinearSolver_t> linear_solver = std::make_shared<OmniLinearSolver_t>();

            // std::shared_ptr<NewtonBase_t> nonlinear_solver = std::make_shared<Newton<Matrix_t>>();
            // nonlinear_solver->set_linear_solver(linear_solver);

            auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix_t, Vector_t, HOMEMADE>>();
            std::shared_ptr<NewtonBase_t> nonlinear_solver = std::make_shared<TrustRegion<Matrix_t>>(subproblem);

            // auto qp_solver = std::make_shared<utopia::MPRGP<PetscMatrix, PetscVector> >();

            in.get("solver", *nonlinear_solver);

            ////////////////////////////////////////////////////////////////////////////////////////////////

            bool zero_initial_guess = input_state.empty();
            Vector_t x;
            if (zero_initial_guess) {
                space.create_vector(x);
                x.set(0.);
            } else {
                x = input_state.data();
            }

            Matrix_t hessian;
            space.create_matrix(hessian);

            Vector_t grad;
            space.create_vector(grad);

            bool ok = true;
            if (zero_initial_guess) {
                material->hessian_and_gradient(x, hessian, grad);

                // We need the negative gradient
                grad = -grad;

                space.apply_constraints(grad);

                ok = linear_solver->solve(hessian, grad, x);

                space.apply_zero_constraints(grad);

                if (debug) print_norms(hessian, grad, x);
            }

            space.apply_constraints(x);

            if (false) {
                // Use this branch for debugging materials
                Vector_t correction;
                space.create_vector(correction);
                correction.set(0.);

                for (int i = 0; i < 4 && ok; ++i) {
                    ok = material->hessian_and_gradient(x, hessian, grad);
                    ok = linear_solver->solve(hessian, grad, correction);

                    x -= correction;

                    if (debug) print_norms(hessian, grad, x);
                }

            } else {
                ok = nonlinear_solver->solve(*material, x);
            }

            space.write(output, x);

            if (!ok) {
                utopia::out() << "Yo!\n";
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_BDDSQP_APP_HPP
