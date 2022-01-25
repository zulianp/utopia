#include "utopia_petsc_Library.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "utopia_petsc_build_ksp.hpp"

#ifdef UTOPIA_WITH_SLEPC
#include "slepcsys.h"
#endif

namespace utopia {
    void PetscLibrary::init(int argc, char *argv[]) {
#ifdef UTOPIA_WITH_PETSC
        static char help[] = "initializing utopia environment through petsc";

        // #ifndef NDEBUG
        // PetscOptionsSetValue(nullptr, "-on_error_attach_debugger", "lldb");
        // PetscOptionsSetValue(nullptr, "-start_in_debugger", "lldb");
        // #else
        PetscOptionsSetValue(nullptr, "-on_error_abort", nullptr);
        // #endif

#ifdef UTOPIA_WITH_SLEPC
        SlepcInitialize(&argc, &argv, (char *)nullptr, help);  // calls PetscInitialize inside
#else
        PetscInitialize(&argc, &argv, (char *)0, help);
#endif  // UTOPIA_WITH_SLEPC

        // is this proper place for doing this ???
        KSPRegister("utopia", KSPCreate_UTOPIA);
#endif
    }

    int PetscLibrary::finalize() {
#ifdef UTOPIA_WITH_SLEPC
        return SlepcFinalize();  // calls PetscFinalize inside
#else
        return PetscFinalize();
#endif  // UTOPIA_WITH_SLEPC
    }

    std::string PetscLibrary::name() const { return "Petsc"; }

    std::string PetscLibrary::version_string() const {
        return std::to_string(PETSC_VERSION_MAJOR) + "." + std::to_string(PETSC_VERSION_MINOR) + "." +
               std::to_string(PETSC_VERSION_SUBMINOR);
    }
}  // namespace utopia
