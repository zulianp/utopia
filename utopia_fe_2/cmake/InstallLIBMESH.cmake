

if(NOT LIBMESH_FOUND)

    #  git clone git://github.com/libMesh/libmesh.git
    include(ExternalProject)

    set(STAGE_DIR            "${CMAKE_BINARY_DIR}/stage")
    set(LIBMESH_URL            git://github.com/libMesh/libmesh.git )
    set(LIBMESH_SOURCE_DIR     ${STAGE_DIR}/libmesh)
    set(LIBMESH_BIN_DIR        ${STAGE_DIR}/libmesh/bin)
    set(LIBMESH_INSTALL_DIR    ${CMAKE_SOURCE_DIR}/external/libmesh)
    set(LIBMESH_MPI_BASE_DIR   $ENV{MPI_DIR})
    set(MAKE_COMMAND "make")

    set(LIBMESH_CONFIG_ARGS "${LIBMESH_CONFIG_ARGS}")

    set(METHODS "opt dbg")

    ExternalProject_Add(
        libmesh
        UPDATE_COMMAND      "" #FIXME
        BUILD_IN_SOURCE 1
        PREFIX              ${STAGE_DIR}
        GIT_REPOSITORY      ${LIBMESH_URL}
	#GIT_TAG		    tags/v1.3.1
        DOWNLOAD_DIR        ${STAGE_DIR}
        INSTALL_DIR         ${LIBMESH_INSTALL_DIR}
        LOG_CONFIGURE       1
        LOG_BUILD           1
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-methods=${METHODS} --enable-silent-rules --enable-unique-id --disable-warnings --disable-maintainer-mode
        BUILD_COMMAND make
        INSTALL_COMMAND     make install
        # COMMAND       ${MAKE_COMMAND}
    )

    set_target_properties(libmesh PROPERTIES EXCLUDE_FROM_ALL TRUE)

    set(LIBMESH_DIR ${LIBMESH_INSTALL_DIR})
    set(LIBMESH_DIR ${LIBMESH_INSTALL_DIR} PARENT_SCOPE)
endif()
