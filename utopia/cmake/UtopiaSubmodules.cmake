if(UTOPIA_ENABLE_ISOLVER)

  find_path(
    ISOLVER_LSOLVE_DIR
    NAMES isolver_lsolve.h
    HINTS ${ISOLVER_DIR}/interfaces/lsolve
          $ENV{ISOLVER_DIR}/isolver/interfaces/lsolve
          ${CMAKE_CURRENT_SOURCE_DIR}/external/isolver/interfaces/lsolve)

  if(NOT ISOLVER_LSOLVE_DIR)
    message(FATAL_ERROR "${ISOLVER_LSOLVE_DIR}")
  endif()

  find_path(
    ISOLVER_NLSOLVE_DIR
    NAMES isolver_function.h
    HINTS ${ISOLVER_DIR}/interfaces/nlsolve
          $ENV{ISOLVER_DIR}/isolver/interfaces/nlsolve
          ${CMAKE_CURRENT_SOURCE_DIR}/external/isolver/interfaces/nlsolve)

  if(NOT ISOLVER_NLSOLVE_DIR)
    message(FATAL_ERROR "${ISOLVER_NLSOLVE_DIR}")
  endif()

  scan_directories(${ISOLVER_LSOLVE_DIR} "." UTOPIA_BUILD_INCLUDES
                   UTOPIA_HEADERS UTOPIA_SOURCES)
  set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES})

  set(UTOPIA_HEADERS ${UTOPIA_HEADERS})
  set(UTOPIA_SOURCES ${UTOPIA_SOURCES})

  scan_directories(${ISOLVER_NLSOLVE_DIR} "." UTOPIA_BUILD_INCLUDES
                   UTOPIA_HEADERS UTOPIA_SOURCES)
  set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES})

  set(UTOPIA_HEADERS ${UTOPIA_HEADERS})
  set(UTOPIA_SOURCES ${UTOPIA_SOURCES})

  add_subdirectory(plug_in_and_out/isolver)

endif()

if(UTOPIA_ENABLE_RAPIDJSON)

endif()

if(UTOPIA_ENABLE_TINY_EXPR)
  set(EXTERNAL_DIR external)

  find_path(
    TINY_EXPR_DIR
    NAMES tinyexpr.h
    HINTS ${EXTERNAL_DIR}/tinyexpr ${TINY_EXPR_DIR} $ENV{TINY_EXPR_DIR}
          ${INSTALL_DIR}/tinyexpr $ENV{INSTALL_DIR}/tinyexpr)

  if(TINY_EXPR_DIR)

    # Add headers and sources to global variables.
    scan_directories(${TINY_EXPR_DIR} "." UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS
                     UTOPIA_SOURCES)
    set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES})

    set(UTOPIA_HEADERS ${UTOPIA_HEADERS})

    set(UTOPIA_SOURCES ${UTOPIA_SOURCES})

    set(UTOPIA_ENABLE_TINY_EXPR ON)
    set(UTOPIA_ADDITIONAL_COMPONENTS ";tinyexpr")
  endif()

endif()

# if(UTOPIA_ENABLE_PARMETIS)

# endif()

# if(UTOPIA_ENABLE_MATRIX_IO)

# endif()
