#include "utopia_petsc_build_ksp.hpp"

#include <utility>
#include "utopia_Instance.hpp"

#undef __FUNCT__
#define __FUNCT__ "KSPSetUp_UTOPIA"
static PetscErrorCode KSPSetUp_UTOPIA(KSP /*ksp*/) {
    PetscFunctionBegin;

    // KSP_UTOPIA         *utopia_Pt = (KSP_UTOPIA*)ksp->data;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSolve_UTOPIA"
static PetscErrorCode KSPSolve_UTOPIA(KSP ksp) {
    PetscFunctionBegin;
    // PetscErrorCode ierr;

    Mat Amat, Pmat;

    PCGetOperators(ksp->pc, &Amat, &Pmat);

    auto *utopia_ls = static_cast<KSP_UTOPIA *>(ksp->data);

    utopia_ls->utopia_set_tolerances(ksp->rtol, ksp->abstol, ksp->divtol, ksp->max_it);
    utopia_ls->utopia_solve_routine(Amat, Pmat, ksp->vec_rhs, ksp->vec_sol);

    // petsc preconditioner combined with utopia solver ...
    utopia_ls->get_convergence_reason(ksp->its, ksp->reason);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetSolveRoutine_UTOPIA"
PetscErrorCode KSPSetSolveRoutine_UTOPIA(
    KSP ksp,
    std::function<void(const Mat &, const Mat &, const Vec &, Vec &)> solve_routine) {
    PetscFunctionBegin;

    auto *utopia_ls = static_cast<KSP_UTOPIA *>(ksp->data);
    utopia_ls->utopia_solve_routine = std::move(solve_routine);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetTolerances_UTOPIA"
PetscErrorCode KSPSetTolerances_UTOPIA(
    KSP ksp,
    std::function<void(const PetscReal &, const PetscReal &, const PetscReal &, const PetscInt &)> set_tolerances) {
    PetscFunctionBegin;

    auto *utopia_ls = static_cast<KSP_UTOPIA *>(ksp->data);
    utopia_ls->utopia_set_tolerances = std::move(set_tolerances);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetGetConvergenceReason_UTOPIA"
PetscErrorCode KSPSetGetConvergenceReason_UTOPIA(
    KSP ksp,
    std::function<void(PetscInt &, KSPConvergedReason &)> convergence_reason) {
    PetscFunctionBegin;

    auto *utopia_ls = static_cast<KSP_UTOPIA *>(ksp->data);
    utopia_ls->get_convergence_reason = std::move(convergence_reason);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPDestroy_UTOPIA"
PetscErrorCode KSPDestroy_UTOPIA(KSP ksp) {
    PetscFunctionBegin;
    // KSP_UTOPIA         *utopia = (KSP_UTOPIA*)ksp->data;

    KSPDestroyDefault(ksp);

    PetscObjectComposeFunction((PetscObject)ksp, "KSPUTOPIASetSolveRoutine_C", NULL);
    PetscObjectComposeFunction((PetscObject)ksp, "KSPSetTolerancesRoutine_C", NULL);
    PetscObjectComposeFunction((PetscObject)ksp, "KSPSetGetConvergenceReason_C", NULL);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPView_UTOPIA"
PetscErrorCode KSPView_UTOPIA(KSP /*ksp*/, PetscViewer viewer) {
    PetscFunctionBegin;

    PetscBool iascii;
    PetscObjectTypeCompare(reinterpret_cast<PetscObject>(viewer), PETSCVIEWERASCII, &iascii);
    // KSP_UTOPIA *utopia_ls = (KSP_UTOPIA*)ksp->data;

    // maybe som future printouts ...
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPSetFromOptions_UTOPIA"
PetscErrorCode KSPSetFromOptions_UTOPIA(PetscOptionItems * /*PetscOptionsObject*/, KSP /*ksp*/) {
    PetscFunctionBegin;
    // PetscErrorCode ierr;

    // not much to add so far...
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "KSPCreate_UTOPIA"
PETSC_EXTERN PetscErrorCode KSPCreate_UTOPIA(KSP ksp) {
    PetscErrorCode ierr;
    KSP_UTOPIA *utopia_solver;

    PetscFunctionBegin;
    ierr = PetscNewLog(ksp, &utopia_solver);
    CHKERRQ(ierr);

    ksp->data = (void *)utopia_solver;

    ksp->ops->setup = KSPSetUp_UTOPIA;
    ksp->ops->solve = KSPSolve_UTOPIA;
    ksp->ops->destroy = KSPDestroy_UTOPIA;
    ksp->ops->view = KSPView_UTOPIA;
    ksp->ops->setfromoptions = KSPSetFromOptions_UTOPIA;

    ksp->ops->buildsolution = KSPBuildSolutionDefault;
    ksp->ops->buildresidual = KSPBuildResidualDefault;

    ierr = PetscObjectComposeFunction((PetscObject)ksp, "KSPUTOPIASetSolveRoutine_C", KSPSetSolveRoutine_UTOPIA);
    CHKERRQ(ierr);
    ierr = PetscObjectComposeFunction((PetscObject)ksp, "KSPSetTolerancesRoutine_C", KSPSetTolerances_UTOPIA);
    CHKERRQ(ierr);
    ierr =
        PetscObjectComposeFunction((PetscObject)ksp, "KSPSetGetConvergenceReason_C", KSPSetGetConvergenceReason_UTOPIA);
    CHKERRQ(ierr);

    m_utopia_warning_once(
        "> FIXME: setting KSPSetSupportedNorm(ksp, KSP_NORM_NATURAL, PC_RIGHT, 1) see:\n"
        "  http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetSupportedNorm.html");

    ierr = KSPSetSupportedNorm(ksp, KSP_NORM_NATURAL, PC_RIGHT, 1);
    CHKERRQ(ierr);

    // Maybe use this?
    // ierr = KSPSetSupportedNorm(ksp, KSP_NORM_UNPRECONDITIONED, 1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

namespace utopia {

    void build_ksp(const std::shared_ptr<utopia::LinearSolver<PetscMatrix, PetscVector>> &lin_solver, KSP &ksp) {
        using namespace utopia;
        using Matrix = utopia::PetscMatrix;
        using Vector = utopia::PetscVector;

        assert(lin_solver);
        assert(ksp);

        //////////////// petsc specializations //////////////////////////

        auto factorization = std::dynamic_pointer_cast<Factorization<Matrix, Vector, PETSC>>(lin_solver);

        if (factorization && factorization->is_proxy()) {
            factorization->strategy().set_ksp_options(ksp);
            return;
        }

        if (!factorization) {
            auto ksp_solver = std::dynamic_pointer_cast<KSPSolver<Matrix, Vector, PETSC>>(lin_solver);
            if (ksp_solver) {
                ksp_solver->set_ksp_options(ksp);
                ksp_solver->attach_preconditioner(ksp);
                return;
            }
        }

        //////////////// All the others //////////////////////////

        // check if our options overwrite this
        KSPSetFromOptions(ksp);
        KSPSetType(ksp, KSPUTOPIA);

        auto iterative_solver = std::dynamic_pointer_cast<utopia::IterativeSolver<Matrix, Vector>>(lin_solver);

        if (iterative_solver) {
            ksp->rtol = iterative_solver->rtol();
            ksp->abstol = iterative_solver->atol();
            ksp->max_it = iterative_solver->max_it();
        }

        std::function<void(const Mat &, const Mat &, const Vec &, Vec &)> solve_routine =
            [lin_solver](const Mat &A, const Mat & /*P*/, const Vec &b, Vec &x) {
                assert(lin_solver);

                Vector b_ut, x_ut;

                // we need to get some better way how to do this
                utopia::convert(x, x_ut);
                utopia::convert(b, b_ut);

                Matrix A_ut;
                A_ut.wrap(const_cast<Mat &>(A));

                lin_solver->solve(A_ut, b_ut, x_ut);

                convert(x_ut, x);
            };

        std::function<void(const PetscReal &, const PetscReal &, const PetscReal &, const PetscInt &)> tol_routine =
            [lin_solver](
                const PetscReal &rtol, const PetscReal &abstol, const PetscReal & /*dtol*/, const PetscInt &maxits) {
                assert(lin_solver);

                if (dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get()) != nullptr) {
                    auto ls = dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get());
                    ls->atol(abstol);
                    ls->rtol(rtol);
                    ls->max_it(maxits);
                }
            };

        std::function<void(PetscInt &, KSPConvergedReason &)> convergence_routine =
            [lin_solver](PetscInt &max_it, KSPConvergedReason &reason) {
                assert(lin_solver);

                if (dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get()) != nullptr) {
                    auto ls = dynamic_cast<utopia::IterativeSolver<Matrix, Vector> *>(lin_solver.get());
                    max_it = ls->get_num_it();

                    switch (ls->get_convergence_reason()) {
                            // sucess
                        case utopia::ConvergenceReason::CONVERGED_FNORM_ABS:
                            reason = KSP_CONVERGED_ATOL;
                            break;

                        case utopia::ConvergenceReason::CONVERGED_FNORM_RELATIVE:
                            reason = KSP_CONVERGED_RTOL;
                            break;

                        case utopia::ConvergenceReason::CONVERGED_SNORM_RELATIVE:
                            reason = KSP_CONVERGED_STEP_LENGTH;
                            break;

                            // fail
                        case utopia::ConvergenceReason::DIVERGED_MAX_IT:
                            reason = KSP_DIVERGED_ITS;
                            break;

                        default:
                            reason = KSP_CONVERGED_RTOL_NORMAL;
                    }
                } else {
                    m_utopia_warning_once("get convergence reason is not properly handled for direct solvers yet...");
                    reason = KSP_CONVERGED_RTOL_NORMAL;

                    // #if UTOPIA_PETSC_VERSION_LESS_THAN(3, 11, 0)
                    //                     reason = KSP_DIVERGED_PCSETUP_FAILED;
                    // #else
                    //                     reason = KSP_DIVERGED_PC_FAILED;
                    // #endif
                }
            };

        KSPSetSolveRoutine_UTOPIA(ksp, solve_routine);
        KSPSetTolerances_UTOPIA(ksp, tol_routine);
        KSPSetGetConvergenceReason_UTOPIA(ksp, convergence_routine);

        {
            // we could hook up PC to precondition utopia solver in future...
            PC pc;
            KSPGetPC(ksp, &pc);
            PCSetType(pc, "none");
        }

        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }

}  // namespace utopia
