# find dependencies
set(PETSC_TEST_RUNS TRUE)
find_package(Petsc)

if(PETSC_FOUND)
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${PETSC_INCLUDES})
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${PETSC_LIBRARIES})

    set(CMAKE_C_COMPILER ${PETSC_COMPILER})

    if(NOT MPI_CXX_COMPILER)
        set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
        message(STATUS "compiler ${MPI_CXX_COMPILER}")
    endif()

    if(MPI_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
    else()
         execute_process(COMMAND mpicxx -v RESULT_VARIABLE MPICXX_FAILED)

        if(MPICXX_FAILED)
            message(STATUS "Using CMAKE compiler, you can define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
        else()
            message(STATUS "-----------------------------------------------")
            message(STATUS "\n[MPI] using mpicxx for compiling c++ files.\nIf you want to use your own compiler define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
            message(STATUS "-----------------------------------------------")
            set(CMAKE_CXX_COMPILER mpicxx)
            set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
        endif()
    endif()

    set(UTOPIA_WITH_PETSC TRUE)
else()
    message(WARNING "[Warning] Petsc not found")
    set(UTOPIA_WITH_PETSC FALSE)
endif()



if(PETSC_FOUND AND UTOPIA_ENABLE_SLEPC)
    find_package(SLEPc)
    if(SLEPC_FOUND)
        list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${SLEPC_INCLUDES})
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${SLEPC_LIBRARIES})
        message(STATUS "Slepc FOUND")
        set(UTOPIA_WITH_SLEPC TRUE)
    else()
        message(WARNING "[Warning] Slepc not found")
        set(UTOPIA_WITH_SLEPC FALSE)
    endif()
endif()


if(SLEPC_FOUND)
    set(UTOPIA_WITH_SLEPC ON PARENT_SCOPE)
    set(UTOPIA_SLEPC TRUE PARENT_SCOPE)
endif()


# set-up module
if(PETSC_FOUND)
    set(UTOPIA_PETSC TRUE PARENT_SCOPE)
    set(UTOPIA_WITH_PETSC ON PARENT_SCOPE)

    list(APPEND PETSC_MODULES
        .
        eval
        solvers
        solvers/linear
        solvers/multilevel
        solvers/nonlinear
        solvers/smoothers
        solvers/nonlinear/constrained
        solvers/nonlinear/constrained/quadratic_programming
        solvers/nonlinear/trust_region
        solvers/nonlinear/line_search
        )

    if(UTOPIA_WITH_SLEPC)
        list(APPEND PETSC_MODULES
            solvers/eigensolvers
        )
    endif()

    if(ENABLE_PASSO_EXTENSIONS)
        list(APPEND PETSC_MODULES
                external
            )
    endif()

    utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${PETSC_MODULES}")
    target_include_directories(utopia SYSTEM PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(utopia PUBLIC ${PETSC_LIBRARIES})
    target_compile_definitions(utopia PUBLIC ${PETSC_DEFINITIONS})

    if(SLEPC_FOUND)
        target_include_directories(utopia SYSTEM PUBLIC ${SLEPC_INCLUDES})
        target_link_libraries(utopia PUBLIC ${SLEPC_LIBRARIES})
        target_compile_definitions(utopia PUBLIC ${SLEPC_LIB_SLEPC})
    endif()

    set(UTOPIA_THIRDPARTY_INCLUDES ${UTOPIA_THIRDPARTY_INCLUDES} PARENT_SCOPE)
    set(UTOPIA_THIRDPARTY_LIBRARIES ${UTOPIA_THIRDPARTY_LIBRARIES} PARENT_SCOPE)
endif()



