#include "utopia_petsc_ConvergedReason.hpp"

#include "utopia_Config.hpp"
#include "utopia_petsc_Base.hpp"

#include <petscksp.h>

namespace utopia {

    bool ksp_convergence_fatal_failure(const int code) {
        switch (code) {
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 11, 0)
            case KSP_DIVERGED_PCSETUP_FAILED:
#else
            case KSP_DIVERGED_PC_FAILED:
#endif
            {
                return true;
            }

            default: { return false; }
        }
    }

    void print_ksp_converged_reason(const int code, std::ostream &os) {
        switch (code) {
            case 0: {
                os << "CONVERGED\n";
                break;
            }

            case KSP_CONVERGED_RTOL_NORMAL: {
                os << "KSP_CONVERGED_RTOL_NORMAL\n";
                break;
            }

            case KSP_CONVERGED_ATOL_NORMAL: {
                os << "KSP_CONVERGED_ATOL_NORMAL\n";
                break;
            }

            case KSP_CONVERGED_RTOL: {
                os << "KSP_CONVERGED_RTOL\n";
                break;
            }

            case KSP_CONVERGED_ATOL: {
                os << "KSP_CONVERGED_ATOL\n";
                break;
            }

            case KSP_CONVERGED_ITS: {
                os << "KSP_CONVERGED_ITS\n";
                break;
            }

            case KSP_CONVERGED_CG_NEG_CURVE: {
                os << "KSP_CONVERGED_CG_NEG_CURVE\n";
                break;
            }

            case KSP_CONVERGED_CG_CONSTRAINED: {
                os << "KSP_CONVERGED_CG_CONSTRAINED\n";
                break;
            }

            case KSP_CONVERGED_STEP_LENGTH: {
                os << "KSP_CONVERGED_STEP_LENGTH\n";
                break;
            }

            case KSP_CONVERGED_HAPPY_BREAKDOWN: {
                os << "KSP_CONVERGED_HAPPY_BREAKDOWN\n";
                break;
            }

                /* diverged */
            case KSP_DIVERGED_NULL: {
                os << "KSP_DIVERGED_NULL\n";
                break;
            }

            case KSP_DIVERGED_ITS: {
                os << "KSP_DIVERGED_ITS\n";
                break;
            }

            case KSP_DIVERGED_DTOL: {
                os << "KSP_DIVERGED_DTOL\n";
                break;
            }

            case KSP_DIVERGED_BREAKDOWN: {
                os << "KSP_DIVERGED_BREAKDOWN\n";
                break;
            }

            case KSP_DIVERGED_BREAKDOWN_BICG: {
                os << "KSP_DIVERGED_BREAKDOWN_BICG\n";
                break;
            }

            case KSP_DIVERGED_NONSYMMETRIC: {
                os << "KSP_DIVERGED_NONSYMMETRIC\n";
                break;
            }

            case KSP_DIVERGED_INDEFINITE_PC: {
                os << "KSP_DIVERGED_INDEFINITE_PC\n";
                break;
            }

            case KSP_DIVERGED_NANORINF: {
                os << "KSP_DIVERGED_NANORINF\n";
                break;
            }

            case KSP_DIVERGED_INDEFINITE_MAT: {
                os << "KSP_DIVERGED_INDEFINITE_MAT\n";
                break;
            }
#if UTOPIA_PETSC_VERSION_LESS_THAN(3, 11, 0)
            case KSP_DIVERGED_PCSETUP_FAILED:
#else
            case KSP_DIVERGED_PC_FAILED:
#endif
            {
                os << "KSP_DIVERGED_PCSETUP_FAILED\n";
                break;
            }

            default: {
                os << "unknown code\n";
                break;
            }
        }
    }
}  // namespace utopia
