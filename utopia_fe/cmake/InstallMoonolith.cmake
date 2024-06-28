if(UTOPIA_INSTALL_MOONOLITH_EXTERNAL)
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

else()

  include(FetchContent)
  message(STATUS "Fetching par_moonolith, due to option UTOPIA_INSTALL_MOONOLITH=ON.")

  set(MOONOLITH_ENABLE_BENCHMARK
      OFF
      CACHE INTERNAL "")

  set(MOONOLITH_ENABLE_TESTING
      OFF
      CACHE INTERNAL "")

  # FIXME This is used to avoid clashes with test_install and other targets
  set(MOONOLITH_ENABLE_SUBMODULE ON CACHE INTERNAL "")

  FetchContent_Declare(
      moonolith
      GIT_REPOSITORY https://bitbucket.org/zulianp/par_moonolith.git
      # GIT_TAG origin/sampler
      GIT_TAG origin/development
  )
  
  FetchContent_MakeAvailable(moonolith)

  add_library(ParMoonolith::par_moonolith ALIAS par_moonolith)

  list(APPEND UTOPIA_FE_SUBMODULES "par_moonolith")
endif()
