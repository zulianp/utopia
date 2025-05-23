if(UTOPIA_ENABLE_SLEPC)

  if(UTOPIA_DEPENDENCIES_DIR)
    set(SLEPC_INSTALL_DIR ${UTOPIA_DEPENDENCIES_DIR}/SLEPC)
  else()
    set(SLEPC_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../external/SLEPC)
  endif()

  set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
  set(SLEPC_URL https://gitlab.com/slepc/slepc)
  set(SLEPC_SOURCE_DIR ${STAGE_DIR}/slepc)
  set(SLEPC_BIN_DIR ${STAGE_DIR}/slepc/bin)

  set(SLEPC_MPI_BASE_DIR $ENV{MPI_DIR})
  set(MAKE_COMMAND "make")

  set(SLEPC_CONFIG_ARGS $ENV{SLEPC_CONFIG_ARGS})

  ExternalProject_Add(
      slepc
      UPDATE_COMMAND "" # FIXME
      BUILD_IN_SOURCE 1
      PREFIX ${STAGE_DIR}
      GIT_REPOSITORY ${SLEPC_URL}
      GIT_TAG main
      DOWNLOAD_DIR ${STAGE_DIR}
      INSTALL_DIR ${SLEPC_INSTALL_DIR}
      LOG_CONFIGURE 1
      LOG_BUILD 1
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
                        ${SLEPC_CONFIG_ARGS}
      BUILD_COMMAND ${MAKE_COMMAND}
      INSTALL_COMMAND make install
      # COMMAND       ${MAKE_COMMAND}
    )

    set_target_properties(slepc PROPERTIES EXCLUDE_FROM_ALL TRUE)
    set(SLEPC_DIR ${SLEPC_INSTALL_DIR})
    set(ENV{SLEPC_DIR} ${SLEPC_INSTALL_DIR})

endif()
