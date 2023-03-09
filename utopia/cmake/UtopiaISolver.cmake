# UtopiaISolver.cmake

option(UTOPIA_ENABLE_ISOLVER "Enable plug-out interface for isolver" ON)

if(UTOPIA_ENABLE_ISOLVER)
    find_path(
        ISOLVER_LSOLVE_DIR
        NAMES isolver_lsolve.h
        HINTS ${ISOLVER_DIR}/interfaces/lsolve
        	  $ENV{ISOLVER_DIR}/isolver/interfaces/lsolve
        	  ${CMAKE_CURRENT_SOURCE_DIR}/external/isolver/interfaces/lsolve
    )

    if(NOT ISOLVER_LSOLVE_DIR)
        message(FATAL_ERROR "${ISOLVER_LSOLVE_DIR}")
    endif()

    utopia_add_external_library(${ISOLVER_LSOLVE_DIR} ".")

    find_path(
        ISOLVER_NLSOLVE_DIR
        NAMES isolver_function.h
        HINTS ${ISOLVER_DIR}/interfaces/nlsolve
              $ENV{ISOLVER_DIR}/isolver/interfaces/nlsolve
              ${CMAKE_CURRENT_SOURCE_DIR}/external/isolver/interfaces/nlsolve
    )

    if(NOT ISOLVER_NLSOLVE_DIR)
        message(FATAL_ERROR "${ISOLVER_NLSOLVE_DIR}")
    endif()

    utopia_add_external_library(${ISOLVER_NLSOLVE_DIR} ".")

endif()

