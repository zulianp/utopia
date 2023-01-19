#include "utopia_fe_base.hpp"

#include "utopia_Field.hpp"
#include "utopia_MPITimeStatistics.hpp"

#include "utopia_ui.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include "utopia_fe_Core.hpp"

#include "utopia_moonolith.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"

#include "utopia_ContactFEFunction.hpp"
#include "utopia_FEFunctionFactory.hpp"
#include "utopia_NewmarkIntegrator.hpp"
// #include "utopia_moonolith_ContactFEFunctionFactory.hpp"

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

        using ContactFEFunction_t = utopia::ContactFEFunction<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using MaterialFactory_t = utopia::kokkos::MaterialFactory<FunctionSpace, FE>;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        using NewtonBase_t = utopia::NewtonBase<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        void run(Input &in) {
            FunctionSpace space;
            Field<FunctionSpace> deformation;
            bool use_state_as_initial_condition = false;

            in.get("use_state_as_initial_condition", use_state_as_initial_condition);

            in.require("space", [&](Input &in) {
                bool read_state = false;
                in.get("read_state", read_state);

                if (read_state) {
                    space.read_with_state(in, deformation);

                    const Scalar_t norm_deformation = norm2(deformation.data());
                    utopia::out() << "norm_deformation: " << norm_deformation << std::endl;

                } else {
                    space.read(in);
                }
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

            ////////////////////////////////////////////////////////////////////////////////////////////////

            std::shared_ptr<LinearSolver_t> linear_solver = std::make_shared<OmniLinearSolver_t>();
            std::shared_ptr<NewtonBase_t> nonlinear_solver = std::make_shared<Newton<Matrix_t>>();

            nonlinear_solver->set_linear_solver(linear_solver);
            in.get("solver", *nonlinear_solver);

            ////////////////////////////////////////////////////////////////////////////////////////////////

            bool zero_initial_guess = deformation.empty();
            Vector_t x;
            if (zero_initial_guess) {
                space.create_vector(x);
            } else {
                x = deformation.data();
            }

            Matrix_t hessian;
            space.create_matrix(hessian);

            Vector_t grad;
            space.create_vector(grad);

            bool ok = true;
            if (zero_initial_guess) {
                material->hessian_and_gradient(x, hessian, grad);

                // We need the negative gradient
                grad *= -1;

                space.apply_constraints(grad);

                ok = linear_solver->solve(hessian, grad, x);
            }

            if (!ok) {
                utopia::out() << "Yo!\n";
            }

            space.apply_constraints(x);
            nonlinear_solver->solve(*material, x);

            space.write(output, x);
        }
    };

}  // namespace utopia

#include "utopia_Main.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
// #include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_stk_intrepid2_Discretization.hpp"
#include "utopia_stk_intrepid2_Material.hpp"

void stk_bddsqp(utopia::Input &in) {
    utopia::BDDSQPApp<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_bddsqp);
