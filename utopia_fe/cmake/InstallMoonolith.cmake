# if(NOT ParMoonolith_FOUND)
  include(ExternalProject)

  if(UTOPIA_DEPENDENCIES_DIR)
    set(MOONOLITH_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/par_moonolith)
  else()
    set(MOONOLITH_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/ParMoonolith)
    message(STATUS "MOONOLITH_INSTALL_DIR=${MOONOLITH_INSTALL_DIR}")
  endif()

  set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
  set(MOONOLITH_URL https://bitbucket.org/zulianp/par_moonolith.git)
  set(MOONOLITH_SOURCE_DIR ${STAGE_DIR}/moonolith)
  set(MOONOLITH_BIN_DIR ${STAGE_DIR}/moonolith/bin)
  set(MOONOLITH_MPI_BASE_DIR $ENV{MPI_DIR})

  list(APPEND 
    MOONOLITH_CMAKE_ARGS
    "-DCMAKE_INSTALL_PREFIX=${MOONOLITH_INSTALL_DIR}"
    "-DMOONOLITH_ENABLE_TESTING=OFF"
    "-DMOONOLITH_ENABLE_BENCHMARK=OFF")

  ExternalProject_Add(
    par_moonolith
    UPDATE_COMMAND "" # FIXME
    PREFIX ${STAGE_DIR}
    GIT_REPOSITORY ${MOONOLITH_URL}
    DOWNLOAD_DIR ${STAGE_DIR}
    INSTALL_DIR ${MOONOLITH_INSTALL_DIR}
    # BINARY_DIR                      ${MOONOLITH_SOURCE_DIR}
    CMAKE_ARGS "${MOONOLITH_CMAKE_ARGS}"
    LOG_CONFIGURE 1
    LOG_BUILD 1
    BUILD_COMMAND ${CMAKE_COMMAND} -E echo "Starting $<CONFIG> build"
    COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config $<CONFIG>
    COMMAND ${CMAKE_COMMAND} -E echo "$<CONFIG> build complete")

  # set_target_properties(par_moonolith PROPERTIES EXCLUDE_FROM_ALL TRUE)

  list(APPEND UTOPIA_FE_TARGET_DEPENDENCIES par_moonolith)

  set(MOONOLITH_DIR ${MOONOLITH_INSTALL_DIR})
# endif()
