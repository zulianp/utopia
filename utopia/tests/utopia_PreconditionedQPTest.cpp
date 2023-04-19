
#include "utopia.hpp"
#ifdef UTOPIA_WITH_PETSC

#include "utopia_Newton.hpp"
#include "utopia_SubCommUnitTest.hpp"
#include "utopia_Testing.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "utopia_petsc_PMPRGP.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class PreconditionedQPTest : public SubCommUnitTest<Vector> {
    public:
        void run() override {
            UTOPIA_RUN_TEST(mprgp);
            UTOPIA_RUN_TEST(pgs);
            UTOPIA_RUN_TEST(preconditioned_mprgp);
        }

    private:
        using Traits = utopia::Traits<Vector>;

        template <class Problem>
        void solve_and_verify(QPSolver<Matrix, Vector> &solver, Problem &problem) const {
            problem.run(this->comm(), 200000, true, solver, true);
        }

        void mprgp() {
            QPSolverTestProblem<Matrix, Vector> fun;
            MPRGP<Matrix, Vector> mprgp;
            solve_and_verify(mprgp, fun);
        }

        void pgs() {
            QPSolverTestProblem<Matrix, Vector> fun;
            ProjectedGaussSeidel<Matrix, Vector> pgs;
            solve_and_verify(pgs, fun);
        }

        void preconditioned_mprgp() {
            // QPSolverTestProblem<Matrix, Vector> fun;
            // MPRGP<Matrix, Vector> mprgp;
            // auto pgs = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            // pgs->max_it(1);
            // pgs->n_local_sweeps(2);
            // // pgs->verbose(true);

            // mprgp.set_preconditioner(pgs);
            // mprgp.rtol(1e-14);
            // mprgp.stol(1e-14);
            // mprgp.atol(1e-14);

            // solve_and_verify(mprgp, fun);


            QPSolverTestProblem<Matrix, Vector> fun;
            PMPRGP<Matrix, Vector> pmprgp;
            pmprgp.rtol(1e-14);
            pmprgp.stol(1e-14);
            pmprgp.atol(1e-14);

            solve_and_verify(pmprgp, fun);
        }
    };

    void psqp_test() {
        const bool verbose = Utopia::instance().verbose();
        // #ifdef UTOPIA_WITH_BLAS
        //         // Serial backend
        //         run_serial_test<PreconditionedQPTest<BlasMatrixd, BlasVectord>>();
        // #endif  // UTOPIA_WITH_BLAS

#ifdef UTOPIA_WITH_PETSC
        run_parallel_test<PreconditionedQPTest<PetscMatrix, PetscVector>>(verbose);
#endif  // UTOPIA_WITH_PETSC

        // #ifdef UTOPIA_WITH_TRILINOS
        //         run_parallel_test<PreconditionedQPTest<TpetraMatrixd, TpetraVectord>>(verbose);
        // #endif  // UTOPIA_WITH_TRILINOS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(psqp_test);

}  // namespace utopia

#endif  // UTOPIA_WITH_PETSC
