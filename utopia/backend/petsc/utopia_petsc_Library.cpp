#include "utopia_petsc_Library.hpp"
#include "utopia_Base.hpp"

#include "petscsys.h"
#include "utopia_petsc_build_ksp.hpp"

#ifdef WITH_SLEPC
#include "slepcsys.h"
#endif

namespace utopia {
    void PetscLibrary::init(int argc, char *argv[]) {
#ifdef UTOPIA_WITH_PETSC
        static char help[] = "initializing utopia environment through petsc";

        PetscOptionsSetValue(nullptr, "-on_error_abort", nullptr);

#ifdef WITH_SLEPC
        SlepcInitialize(&argc, &argv, (char *)nullptr, help);  // calls PetscInitialize inside
#else
        PetscInitialize(&argc, &argv, (char *)0, help);
#endif  // WITH_SLEPC

        // is this proper place for doing this ???
        KSPRegister("utopia", KSPCreate_UTOPIA);
#endif
    }

    int PetscLibrary::finalize() {
#ifdef WITH_SLEPC
        return SlepcFinalize();  // calls PetscFinalize inside
#else
        return PetscFinalize();
#endif  // WITH_SLEPC
    }

    std::string PetscLibrary::name() const { return "Petsc"; }

    std::string PetscLibrary::version_string() const {
        return std::to_string(PETSC_VERSION_MAJOR) + "." + std::to_string(PETSC_VERSION_MINOR) + "." +
               std::to_string(PETSC_VERSION_SUBMINOR);
    }
}  // namespace utopia
