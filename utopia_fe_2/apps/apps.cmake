

set(UTOPIA_FE_APPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/apps)

list(APPEND APPS_MODULES
    .
    )

if(TARGET utopia_libmesh)
    list(APPEND APPS_MODULES
        libmesh
        )


endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(UTOPIA_FE_APPS_DIR APPS_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(utopia_fe_exec PRIVATE ${LOCAL_SOURCES})
# utopia_link_default_targets(utopia_fe_exec)


target_include_directories(utopia_fe_exec PRIVATE ${UTOPIA_FE_APPS_DIR})
target_include_directories(utopia_fe_exec PRIVATE .)
target_include_directories(utopia_fe_exec PRIVATE ${APPS_MODULES})

if(Gperftools_FOUND)
    target_link_libraries(utopia_fe_exec gperftools::profiler)
endif()

if(TRY_WITH_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(WITH_EIGEN_3 ON)
        set(WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia_fe_exec PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_FE_DEV_FLAGS}")


