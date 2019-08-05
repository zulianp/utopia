#ifndef UTOPIA_PETSC_BASE_HPP
#define UTOPIA_PETSC_BASE_HPP

#include "utopia_Base.hpp"
#include "utopia_Config.hpp"
#include "petscversion.h"

#ifdef WITH_PETSC
    #define UTOPIA_PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
      ((PETSC_VERSION_MAJOR < (major) ||                   \
        (PETSC_VERSION_MAJOR == (major) && ((PETSC_VERSION_MINOR < (minor)) || \
                                                             (PETSC_VERSION_MINOR == (minor) && \
                                                              PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

#define UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(major,minor,subminor)                   \
      ((PETSC_VERSION_MAJOR > (major) ||                   \
        (PETSC_VERSION_MAJOR == (major) && ((PETSC_VERSION_MINOR > (minor)) || \
                                                             (PETSC_VERSION_MINOR == (minor) && \
                                                              PETSC_VERSION_SUBMINOR >= (subminor))))) ? 1 : 0)
#else // WITH_PETSC
    #define UTOPIA_PETSC_VERSION_LESS_THAN(major,minor,subminor) 1
    #define UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(major,minor,subminor) 0
#endif // WITH_PETSC

#endif //UTOPIA_PETSC_BASE_HPP
