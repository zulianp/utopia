#ifndef UTOPIA_PETSC_BASE_HPP
#define UTOPIA_PETSC_BASE_HPP

#include "petscversion.h"
#include "utopia_Base.hpp"
#include "utopia_Config.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#define UTOPIA_PETSC_VERSION_LESS_THAN(major, minor, subminor)                                                        \
    ((PETSC_VERSION_MAJOR < (major) ||                                                                                \
      (PETSC_VERSION_MAJOR == (major) &&                                                                              \
       ((PETSC_VERSION_MINOR < (minor)) || (PETSC_VERSION_MINOR == (minor) && PETSC_VERSION_SUBMINOR < (subminor))))) \
         ? 1                                                                                                          \
         : 0)

#define UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(major, minor, subminor)                                                \
    ((PETSC_VERSION_MAJOR > (major) ||                                                                                 \
      (PETSC_VERSION_MAJOR == (major) &&                                                                               \
       ((PETSC_VERSION_MINOR > (minor)) || (PETSC_VERSION_MINOR == (minor) && PETSC_VERSION_SUBMINOR >= (subminor))))) \
         ? 1                                                                                                           \
         : 0)
#else  // UTOPIA_ENABLE_PETSC
#define UTOPIA_PETSC_VERSION_LESS_THAN(major, minor, subminor) 1
#define UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(major, minor, subminor) 0
#endif  // UTOPIA_ENABLE_PETSC

#endif  // UTOPIA_PETSC_BASE_HPP
