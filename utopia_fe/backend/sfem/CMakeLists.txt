list(APPEND SFEM_MODULES .)

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${SFEM_MODULES}"
                   UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

  set(UTOPIA_FE_BUILD_INCLUDES
      ${UTOPIA_FE_BUILD_INCLUDES}
      PARENT_SCOPE)

  set(UTOPIA_FE_HEADERS
      ${UTOPIA_FE_HEADERS}
      PARENT_SCOPE)

  set(UTOPIA_FE_SOURCES
      ${UTOPIA_FE_SOURCES}
      PARENT_SCOPE)

find_path(
    SFEM_DIR
    NAMES  base/sfem_base.h
    HINTS ${SFEM_DIR} $ENV{SFEM_DIR}/sfem ${CMAKE_SOURCE_DIR}/external/sfem
)

if(SFEM_DIR)
    list(APPEND SFEM_SUBMODULES
        .
        algebra
        base
        baseline
        eigen
        matrix
        mesh
        operators
        operators
        operators/cvfem/
        operators/kernels/
        operators/tet10/
        operators/tet4/
        operators/macro_tet4/
        operators/cvfem_tet4
        operators/tri3/
        operators/macro_tri3/
        operators/tri6/
        operators/trishell3
        operators/beam2/
        perf
        pizzastack
        resampling
        solver
        operators/cvfem_quad4
        operators/cvfem_tet4
        operators/cvfem_tri3
        operators/hierarchical
        operators/kernels
        operators/macro_tet4
        operators/macro_tri3
        operators/navier_stokes
        operators/phase_field_for_fracture
        operators/tet10
        operators/tet4
        operators/tri3
        operators/tri6
        operators/trishell3
    )

    scan_directories(${SFEM_DIR} "${SFEM_SUBMODULES}"
                     UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

    set(UTOPIA_FE_BUILD_INCLUDES
        ${UTOPIA_FE_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_FE_HEADERS
        ${UTOPIA_FE_HEADERS}
        PARENT_SCOPE)

    set(UTOPIA_FE_SOURCES
        ${UTOPIA_FE_SOURCES}
        PARENT_SCOPE)


    set(SFEM_FOUND TRUE)
    set(PARENT_SCOPE SFEM_FOUND TRUE)

    # FIXME Better integration with SFEM in the next versions!

    # utopiafe_add_module(utopia_dep_sfem ${SFEM_DIR} "${SFEM_SUBMODULES}")

    # file(GLOB SFEM_PYTHON_SCRIPTS_1 "${SFEM_DIR}/python/mesh/*.py")
    # file(GLOB SFEM_PYTHON_SCRIPTS_2 "${SFEM_DIR}/python/sdf/*.py")
    # file(GLOB SFEM_PYTHON_SCRIPTS_2 "${SFEM_DIR}/python/grid/*.py")

    # list(APPEND SFEM_PYTHON_SCRIPTS ${SFEM_PYTHON_SCRIPTS_1})
    # list(APPEND SFEM_PYTHON_SCRIPTS ${SFEM_PYTHON_SCRIPTS_2})

    # install(FILES ${SFEM_PYTHON_SCRIPTS} 
    #     PERMISSIONS OWNER_READ OWNER_EXECUTE 
    #     DESTINATION scripts)

    # message(STATUS "SFEM_PYTHON_SCRIPTS=${SFEM_PYTHON_SCRIPTS}")

    # file(GLOB SFEM_EXCUTABLES "${SFEM_DIR}/drivers/*.c")

    # set(SFEM_EXCLUDED_EXEC "")
    # list(APPEND SFEM_EXCLUDED_EXEC "grad_and_project;heat_equation;integrate_divergence;partition_mesh_based_on_operator;taylor_hood_navier_stokes")

    # # message(STATUS "SFEM_EXCLUDED_EXEC=${SFEM_EXCLUDED_EXEC}")

    # # Add sfem executables as optional functionalities
    # foreach(SFEM_EXEC ${SFEM_EXCUTABLES})
    #     get_filename_component(SFEM_EXEC_NAME ${SFEM_EXEC} NAME_WE)
    #     list(FIND SFEM_EXCLUDED_EXEC "${SFEM_EXEC_NAME}" LIST_INDEX)
        
    #     if (${LIST_INDEX} GREATER -1)
    #     else()
    #         # message(STATUS "SFEM_EXEC_NAME=${SFEM_EXEC_NAME}")
    #         # add_executable(${SFEM_EXEC_NAME} EXCLUDE_FROM_ALL ${SFEM_EXEC})
    #         add_executable(${SFEM_EXEC_NAME} ${SFEM_EXEC})
    #         target_link_libraries(${SFEM_EXEC_NAME} PRIVATE utopia_dep_sfem)

    #         install(
    #             TARGETS ${SFEM_EXEC_NAME}
    #             RUNTIME DESTINATION bin
    #             LIBRARY DESTINATION lib
    #             ARCHIVE DESTINATION lib)
    #     endif()
    # endforeach()

    # target_link_libraries(utopia_dep_sfem PUBLIC Utopia::utopia)

    # utopiafe_add_module(utopia_sfem ${CMAKE_CURRENT_SOURCE_DIR} ".")
    # target_link_libraries(utopia_sfem PUBLIC utopia_dep_sfem utopia_fe_core)

    # set(UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES} PARENT_SCOPE)
else()
    message(FATAL_ERROR "Could not find SFEM!")
endif()
