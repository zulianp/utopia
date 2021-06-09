#ifndef UTOPIA_OBSTACLE_APP
#define UTOPIA_OBSTACLE_APP

#include "utopia_fe_base.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_Field.hpp"
#include "utopia_ILU.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_MonotoneAlgebraicMultigrid.hpp"

#include "utopia_MonotoneSemiGeometricMultigrid.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_PatchSmoother.hpp"
#include "utopia_RASPatchSmoother.hpp"

#include "utopia_ui.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include "utopia_fe_Core.hpp"

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

        // Use specialized compoenents for function space
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;

        // Use algorithms from utopia algebra
        using QPSolver_t = utopia::QPSolver<Matrix_t, Vector_t>;
        using SemismoothNewton_t = utopia::SemismoothNewton<Matrix_t, Vector_t>;
        using Factorization_t = utopia::Factorization<Matrix_t, Vector_t>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using IterativeSolver_t = utopia::IterativeSolver<Matrix_t, Vector_t>;
        using AlgebraicMultigrid_t = utopia::AlgebraicMultigrid<Matrix_t, Vector_t>;
        using ProjectedGaussSeidel_t = utopia::ProjectedGaussSeidel<Matrix_t, Vector_t>;
        using KSPSolver_t = utopia::KSPSolver<Matrix_t, Vector_t>;
        using MonotoneAlgebraicMultigrid_t = utopia::MonotoneAlgebraicMultigrid<Matrix_t, Vector_t>;
        using MonotoneSemiGeometricMultigrid_t = utopia::MonotoneSemiGeometricMultigrid<FunctionSpace>;
        using SemiGeometricMultigrid_t = utopia::SemiGeometricMultigridNew<FunctionSpace>;
        using ObstacleFEFunction_t = utopia::ObstacleFEFunction<FunctionSpace>;
        using BoxConstrainedFEFunctionSolver_t = utopia::BoxConstrainedFEFunctionSolver<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;

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

            auto fun = std::make_shared<FEModelFunction_t>(make_ref(space));
            auto obs_fun = std::make_shared<ObstacleFEFunction_t>(fun);

            obs_fun->read(in);

            BoxConstrainedFEFunctionSolver_t solver;
            in.get("solver", solver);

            Vector_t x;
            space.create_vector(x);

            if (!solver.solve(*obs_fun, x)) {
                space.comm().root_print("ObstacleApp[Warning] solver failed to converge!");
            }

            obs_fun->report_solution(x);
        }
    };

}  // namespace utopia

#endif
