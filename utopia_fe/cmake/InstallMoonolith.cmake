if(NOT ParMoonolith_FOUND)
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
    set(MAKE_COMMAND "make")

    set(MOONOLITH_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${MOONOLITH_INSTALL_DIR}")

    set(METHODS "opt dbg")

    ExternalProject_Add(
        par_moonolith
        UPDATE_COMMAND "" # FIXME
        BUILD_IN_SOURCE 1
        PREFIX ${STAGE_DIR}
        GIT_REPOSITORY ${MOONOLITH_URL}
        # GIT_TAG             tags/v1.3.1
        DOWNLOAD_DIR ${STAGE_DIR}
        INSTALL_DIR ${MOONOLITH_INSTALL_DIR}
        CMAKE_ARGS "${MOONOLITH_CMAKE_ARGS}"
        LOG_CONFIGURE 1
        LOG_BUILD 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        # COMMAND       ${MAKE_COMMAND}
    )

    set_target_properties(par_moonolith PROPERTIES EXCLUDE_FROM_ALL TRUE)

    set(MOONOLITH_DIR ${MOONOLITH_INSTALL_DIR})
endif()