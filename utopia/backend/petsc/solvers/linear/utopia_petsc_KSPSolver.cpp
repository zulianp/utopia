#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Types.hpp"

#include <petscksp.h>
#include <petscpc.h>

#include "utopia_petsc_KSPSolver_impl.hpp"

namespace utopia {

    /**
     * @brief      Interface between petsc shell preconditioner and utopia solvers/preconditioners.
     */
    PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y) {
        Preconditioner<PetscVector> *shell;
        PCShellGetContext(pc, reinterpret_cast<void **>(&shell));

        // TODO(zulianp): ref would be nicer here
        PetscVector x_u, y_u;
        convert(x, x_u);
        convert(y, y_u);

        shell->apply(x_u, y_u);
        convert(y_u, y);

        return 0;
    }

    PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void *ctx) {
        auto *logger = static_cast<KSPLog *>(ctx);

        Vec x;
        PetscReal conv_rate = 0.0;

        KSPBuildSolution(ksp, nullptr, &x);

        MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(x));

        if (!logger->initialized()) {
            logger->init_from(x);
        }

        if (it > 1) {
            Vec nom, denom;

            VecDuplicate(x, &denom);
            VecDuplicate(x, &nom);

            VecCopy(logger->x_k_1, nom);
            VecCopy(x, denom);

            VecAYPX(nom, -1, x);
            VecAYPX(denom, -1, logger->x_k_2);

            PetscReal next, previous;
            VecNorm(nom, NORM_2, &next);
            VecNorm(denom, NORM_2, &previous);

            conv_rate = next / previous;

            VecDestroy(&nom);
            VecDestroy(&denom);
        }

        if (it > 0) {
            VecCopy(logger->x_k_1, logger->x_k_2);
        }

        // if(it > 0)
        VecCopy(x, logger->x_k_1);

        PetscBool compute_cond_number;
        KSPGetComputeSingularValues(ksp, &compute_cond_number);

        if (compute_cond_number != 0u) {
            if (it == 0) {
                PetscPrintf(comm, "it           ||r||                   rho                  cond. number \n");
            }

            PetscReal emax, emin;
            KSPComputeExtremeSingularValues(ksp, &emax, &emin);
            PetscPrintf(comm,
                        "%D     %14.12e         %14.12e        %14.12e \n",
                        it,
                        rnorm,
                        conv_rate,
                        std::abs(emax) / std::abs(emin));
        } else {
            if (it == 0) {
                PetscPrintf(comm, "it           ||r||                   rho      \n");
            }

            PetscPrintf(comm, "%D     %14.12e         %14.12e \n", it, rnorm, conv_rate);
        }

        return 0;
    }

    std::string converged_str(KSPConvergedReason reason) {
        switch (reason) {
            /* converged */
            case KSP_CONVERGED_RTOL_NORMAL: {
                return "KSP_CONVERGED_RTOL_NORMAL";
            }
            case KSP_CONVERGED_ATOL_NORMAL: {
                return "KSP_CONVERGED_ATOL_NORMAL";
            }
            case KSP_CONVERGED_RTOL: {
                return "KSP_CONVERGED_RTOL";
            }
            case KSP_CONVERGED_ATOL: {
                return "KSP_CONVERGED_ATOL";
            }
            case KSP_CONVERGED_ITS: {
                return "KSP_CONVERGED_ITS";
            }
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 19, 0)
            case KSP_CONVERGED_CG_NEG_CURVE: {
                return "KSP_CONVERGED_CG_NEG_CURVE";
            }
            // duplicate of KSP_CONVERGED_STEP_LENGTH (since version 3.19)
            case KSP_CONVERGED_CG_CONSTRAINED: {
                return "KSP_CONVERGED_CG_CONSTRAINED";
            }
#endif
            case KSP_CONVERGED_STEP_LENGTH: {
                return "KSP_CONVERGED_STEP_LENGTH";
            }
            case KSP_CONVERGED_HAPPY_BREAKDOWN: {
                return "KSP_CONVERGED_HAPPY_BREAKDOWN";
            }
            /* diverged */
            case KSP_DIVERGED_NULL: {
                return "KSP_DIVERGED_NULL";
            }
            case KSP_DIVERGED_ITS: {
                return "KSP_DIVERGED_ITS";
            }
            case KSP_DIVERGED_DTOL: {
                return "KSP_DIVERGED_DTOL";
            }
            case KSP_DIVERGED_BREAKDOWN: {
                return "KSP_DIVERGED_BREAKDOWN";
            }
            case KSP_DIVERGED_BREAKDOWN_BICG: {
                return "KSP_DIVERGED_BREAKDOWN_BICG";
            }
            case KSP_DIVERGED_NONSYMMETRIC: {
                return "KSP_DIVERGED_NONSYMMETRIC";
            }
            case KSP_DIVERGED_INDEFINITE_PC: {
                return "KSP_DIVERGED_INDEFINITE_PC";
            }
            case KSP_DIVERGED_NANORINF: {
                return "KSP_DIVERGED_NANORINF";
            }
            case KSP_DIVERGED_INDEFINITE_MAT: {
                return "KSP_DIVERGED_INDEFINITE_MAT";
            }
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 11, 0)
            case KSP_DIVERGED_PCSETUP_FAILED: {
                return "KSP_DIVERGED_PCSETUP_FAILED";
            }
#else
            case KSP_DIVERGED_PC_FAILED: {
                return "KSP_DIVERGED_PCSETUP_FAILED";
            }
#endif
            case KSP_CONVERGED_ITERATING: {
                return "KSP_CONVERGED_ITERATING";
            }

            default: {
                return "unknown code";
            }
        }
    }

    // explicit
    template class KSPSolver<PetscMatrix, PetscVector, PETSC>;
    // FIXME
    // template class KSPSolver<PetscMatrix, PetscVector, PETSC>;
}  // namespace utopia
