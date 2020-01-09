list(APPEND POLYMORPHIC_MODULES
    .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR POLYMORPHIC_MODULES LOCAL_HEADERS LOCAL_SOURCES)
add_library(utopia_polymorphic STATIC ${LOCAL_SOURCES})
utopia_link_default_targets(utopia_polymorphic)

target_link_libraries(utopia_polymorphic utopia_core utopia_utils utopia_solvers utopia_differentiable utopia_ui utopia_interoperability)
target_include_directories(utopia_polymorphic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(utopia_polymorphic PUBLIC ..)
target_include_directories(utopia_polymorphic PUBLIC ../backend ${POLYMORPHIC_MODULES})

if(TRY_WITH_EIGEN_3)
    find_package(Eigen3)
    if(EIGEN3_FOUND)
        set(WITH_EIGEN_3 ON)
        set(WITH_EIGEN_3 ON PARENT_SCOPE)
        target_include_directories(utopia_polymorphic PRIVATE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(TARGET utopia_opencl)
    target_link_libraries(utopia_polymorphic utopia_opencl)
endif()

if(TARGET utopia_petsc)
    target_link_libraries(utopia_polymorphic utopia_petsc)
endif()

if(TARGET utopia_blas)
    target_link_libraries(utopia_polymorphic utopia_blas)
endif()

if(TARGET utopia_cuda_cxx)
    target_link_libraries(utopia_polymorphic utopia_cuda_cxx)
endif()

if(TARGET utopia_trilinos)
    target_link_libraries(utopia_polymorphic utopia_trilinos)
endif()

if(TARGET utopia_m3elinsol)
    target_link_libraries(utopia_polymorphic utopia_m3elinsol)
endif()

install(TARGETS utopia_polymorphic
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_polymorphic" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")
