#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

#include "utopia_BarrierMultigrid.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_WITH_PETSC

namespace utopia {

    template <class Matrix, class Vector>
    class BarrierMultigridTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void run() {
            auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();
            auto preconditioner = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            linear_solver->set_preconditioner(preconditioner);

            BarrierMultigrid<Matrix, Vector> mg(linear_solver);
        }
    };

    static void barrier_mg() {
#ifdef UTOPIA_WITH_PETSC
        BarrierMultigridTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
        BarrierMultigridTest<TpetraMatrixd, TpetraVectord>().run();
#endif  // UTOPIA_WITH_TRILINOS

#ifdef UTOPIA_WITH_BLAS
        BarrierMultigridTest<BlasMatrixd, BlasVectord>()
            .run();  // TODO(zulianp): : because blas is missing min operation ....
#endif               // UTOPIA_WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(barrier_mg);
}  // namespace utopia