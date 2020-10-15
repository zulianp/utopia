if(UTOPIA_PETSC)

    set(UTOPIA_APPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/apps)

    list(
        APPEND
        APPS_MODULES
        .
        data_structures
        gpu_assembler
        cpu_assembler
        cpu_assembler/dmplex
        PF_frac
        fe)

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(${UTOPIA_APPS_DIR} "${APPS_MODULES}" LOCAL_HEADERS
                       LOCAL_SOURCES)
    target_sources(utopia_exec PRIVATE ${LOCAL_SOURCES})
    utopia_link_default_targets(utopia_exec)

    target_include_directories(utopia_exec PRIVATE ${UTOPIA_APPS_DIR})
    target_include_directories(utopia_exec PRIVATE .)
    foreach(MODULE ${APPS_MODULES})
        target_include_directories(utopia_exec
                                   PRIVATE ${UTOPIA_APPS_DIR}/${MODULE})
    endforeach()

    if(Gperftools_FOUND)
        target_link_libraries(utopia_exec PUBLIC gperftools::profiler)
    endif()

    if(UTOPIA_ENABLE_EIGEN_3)
        find_package(Eigen3)
        if(EIGEN3_FOUND)
            set(UTOPIA_WITH_EIGEN_3 ON)
            target_include_directories(utopia_exec
                                       PRIVATE ${EIGEN3_INCLUDE_DIR})

        endif()
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")

endif()
