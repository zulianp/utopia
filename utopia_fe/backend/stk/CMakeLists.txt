find_package(Trilinos NO_MODULE NO_CMAKE_ENVIRONMENT_PATH
             NO_SYSTEM_ENVIRONMENT_PATH HINTS $ENV{TRILINOS_DIR})

if(Trilinos_FOUND)
    utopiafe_add_module(utopia_stk ${CMAKE_CURRENT_SOURCE_DIR} ".")

    target_include_directories(utopia_stk PUBLIC ${Trilinos_INCLUDE_DIRS})
    target_include_directories(utopia_stk PUBLIC ${Trilinos_TPL_INCLUDE_DIRS})

    target_link_libraries(utopia_stk PUBLIC Utopia::utopia utopia_fe_core)
    target_link_libraries(utopia_stk INTERFACE ${Trilinos_TPL_LIBRARIES})
    target_link_libraries(utopia_stk INTERFACE ${Trilinos_LIBRARIES})
    target_link_directories(utopia_stk INTERFACE ${Trilinos_TPL_LIBRARY_DIRS})

    set(Trilinos_LIBRARIES
        ${Trilinos_LIBRARIES}
        PARENT_SCOPE)

    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS} ${EXTRA_CXX_FLAGS}")

    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS} ${EXTRA_Fortran_FLAGS}"
    )

    set(CMAKE_CXX_FLAGS
        ${CMAKE_CXX_FLAGS}
        PARENT_SCOPE)

    set(CMAKE_Fortran_FLAGS
        ${CMAKE_Fortran_FLAGS}
        PARENT_SCOPE)

    message(STATUS "[Status] found Trilinos ${Trilinos_INCLUDE_DIRS}")

    set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
    set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

    set(CMAKE_CXX_COMPILER
        ${CMAKE_CXX_COMPILER}
        PARENT_SCOPE)

    set(CMAKE_C_COMPILER
        ${CMAKE_C_COMPILER}
        PARENT_SCOPE)

    set(CMAKE_Fortran_COMPILER
        ${CMAKE_Fortran_COMPILER}
        PARENT_SCOPE)

    set(UTOPIA_WITH_STK TRUE)
    set(UTOPIA_WITH_STK
        TRUE
        PARENT_SCOPE)

else()
    message(WARNING "[Warning] did not find Trilinos")
endif()
