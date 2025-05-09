#ifndef UTOPIA_OBSTACLE_APP
#define UTOPIA_OBSTACLE_APP

#include "utopia_fe_base.hpp"

#include "utopia_Field.hpp"
#include "utopia_MPITimeStatistics.hpp"

#include "utopia_ui.hpp"

#ifdef UTOPIA_ENABLE_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_ENABLE_BLAS

#include "utopia_fe_Core.hpp"

#include "utopia_FEFunctionFactory.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_ObstacleFEFunction.hpp"
#include "utopia_moonolith_ObstacleFEFunctionFactory.hpp"

#include "utopia_MeshTransform.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class ObstacleApp {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        using ObstacleFEFunction_t = utopia::ObstacleFEFunction<FunctionSpace>;
        using BoxConstrainedFEFunctionSolver_t = utopia::BoxConstrainedFEFunctionSolver<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;

        void run(Input &in) {
            FunctionSpace space;
            Vector_t x;
            std::shared_ptr<FEFunctionInterface_t> fun;
            std::shared_ptr<ObstacleFEFunction_t> obs_fun;
            bool enable_restart = false;
            in.get("enable_restart", enable_restart);

            if (enable_restart) {
                IO<FunctionSpace> input_db(space);
                input_db.set_import_all_data(true);
                in.require("space", [&](Input &node) { input_db.open_input(node); });

                Scalar_t t = -1;
                in.get("restart_time", t);

                if (t < 0) {
                    input_db.load_last_time_step();
                    t = input_db.max_time();
                } else {
                    input_db.load_time_step(t);
                }

                fun = std::make_shared<FEModelFunction_t>(make_ref(space));
                auto newmark = std::make_shared<NewmarkIntegrator_t>(fun);
                fun = newmark;

                obs_fun = std::make_shared<ObstacleFEFunction_t>(fun);
                obs_fun->read(in);
                newmark->time()->restart(t);
                obs_fun->setup_IVP(input_db);

                // Copy the solution that was loaded from disk
                x = newmark->solution();
            } else {
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

                in.get("displace", [&](Input &in) {
                    MeshTransform<FunctionSpace> transform;
                    transform.read(in);
                    transform.generate_displacement_field(space, deformation);

                    const Scalar_t norm_deformation = norm2(deformation.data());

                    if (mpi_world_rank() == 0) {
                        utopia::out() << "displace (norm_deformation): " << norm_deformation << std::endl;
                    }
                });

                if (space.empty()) {
                    Utopia::Abort("Space is empty!");
                    return;
                }

                fun = std::make_shared<FEModelFunction_t>(make_ref(space));

                std::string integrator;
                in.get("integrator", integrator);

                if (integrator.empty()) {
                    bool dynamic = false;
                    in.get("dynamic", dynamic);

                    if (dynamic) {
                        fun = std::make_shared<NewmarkIntegrator_t>(fun);
                    }

                } else {
                    fun = ObstacleFEFunctionFactory<FunctionSpace>::make_time_integrator(fun, integrator);
                }

                obs_fun = std::make_shared<ObstacleFEFunction_t>(fun);
                obs_fun->read(in);

                if (deformation.empty()) {
                    space.create_vector(x);
                } else {
                    x = deformation.data();
                }

                fun->setup_IVP(x);

                if (use_state_as_initial_condition) {
                    utopia::out() << "Using state as initial condition!\n";
                    assert(!deformation.empty());
                    if (deformation.empty() || !fun->set_initial_condition(deformation.data())) {
                        Utopia::Abort("Called set_initial_condition on function that does not support it!");
                    }
                }
            }

            // Set up obstacle and solve
            BoxConstrainedFEFunctionSolver_t solver;
            in.get("solver", solver);

            bool skip_solve = false;
            in.get("skip_solve", skip_solve);

            if (skip_solve) {
                fun->update_IVP(x);
                obs_fun->report_solution(x);
                return;
            }

            do {
                if (fun->is_time_dependent()) {
                    const Scalar_t t = fun->time()->get();
                    status("time: " + std::to_string(t));
                }

                if (!solver.solve(*obs_fun, x)) {
                    space.comm().root_print("ObstacleApp[Warning] solver failed to converge!");
                }

                fun->update_IVP(x);
                obs_fun->report_solution(x);

            } while (!fun->is_IVP_solved());
        }

        static void status(const std::string &message) {
            if (!mpi_world_rank()) {
                utopia::out() << message << "\n";
            }
        }
    };

}  // namespace utopia

#endif
