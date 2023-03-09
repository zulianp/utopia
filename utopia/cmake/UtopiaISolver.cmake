# UtopiaISolver.cmake

option(UTOPIA_ENABLE_ISOLVER "Enable plug-out interface for isolver" ON)

if(UTOPIA_ENABLE_ISOLVER)
    find_path(
        ISOLVER_DIR
        NAMES isolver_lsolve.h
        HINTS ${ISOLVER_DIR}/interfaces/lsolve
        	  $ENV{ISOLVER_DIR}/isolver/interfaces/lsolve
        	  ${CMAKE_CURRENT_SOURCE_DIR}/external/isolver/interfaces/lsolve
    )

    if(NOT ISOLVER_DIR)
        message(FATAL_ERROR "${ISOLVER_DIR}")
    endif()

    message("ISOLVER_DIR=${ISOLVER_DIR}")
    
    if(ISOLVER_DIR)
        utopia_add_external_library(${ISOLVER_DIR} ".")
    endif()
endif()

