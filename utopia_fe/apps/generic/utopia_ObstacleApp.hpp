#ifndef UTOPIA_OBSTACLE_APP
#define UTOPIA_OBSTACLE_APP

#include "utopia_fe_base.hpp"

#include "utopia_Field.hpp"
#include "utopia_MPITimeStatistics.hpp"

#include "utopia_ui.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include "utopia_fe_Core.hpp"

#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_ObstacleFEFunction.hpp"

#include <memory>

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
            Field<FunctionSpace> deformation;

            in.get("space", [&](Input &in) {
                bool read_state = false;
                in.get("read_state", read_state);

                if (read_state) {
                    space.read_with_state(in, deformation);

                    const Scalar_t norm_deformation = norm2(deformation.data());
                    std::cout << "norm_deformation: " << norm_deformation << std::endl;

                } else {
                    space.read(in);
                }
            });

            if (space.empty()) {
                return;
            }

            std::shared_ptr<FEFunctionInterface_t> fun = std::make_shared<FEModelFunction_t>(make_ref(space));

            bool dynamic = false;
            in.get("dynamic", dynamic);

            if (dynamic) {
                fun = std::make_shared<NewmarkIntegrator_t>(fun);
            }

            auto obs_fun = std::make_shared<ObstacleFEFunction_t>(fun);
            obs_fun->read(in);

            BoxConstrainedFEFunctionSolver_t solver;
            in.get("solver", solver);

            Vector_t x;
            space.create_vector(x);

            fun->setup_IVP(x);

            do {
                if (!solver.solve(*obs_fun, x)) {
                    space.comm().root_print("ObstacleApp[Warning] solver failed to converge!");
                }

                fun->update_IVP(x);
                obs_fun->report_solution(x);

            } while (!fun->is_IVP_solved());
        }
    };

}  // namespace utopia

#endif
