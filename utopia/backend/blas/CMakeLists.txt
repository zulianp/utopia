# find dependencies
if(PITZ_DORA)
    find_package(OpenBLAS)
    if(OPEN_BLAS_FOUND)
        list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${BLAS_INCLUDE_DIR})
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${BLAS_LIBRARIES})
        list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
        set(WITH_OPEN_BLAS ON)
    else()
        message(WARNING "[Warning] blas not found")
        set(WITH_BLAS OFF)
    endif()
else()
    SET(OPEN_BLAS_FOUND FALSE)
    find_package(BLAS)
    if(BLAS_FOUND)
        list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${BLAS_INCLUDE_DIR})
        list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${BLAS_LIBRARIES})
        list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
        set(WITH_BLAS ON)
    else()
        SET(BLAS_FOUND FALSE)
    endif()
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${LAPACK_INCLUDE_DIR})
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${LAPACK_DEFINITIONS})
    set(WITH_LAPACK ON)
else()
    message(WARNING "[Warning] lapack not found")
    set(WITH_LAPACK OFF)
endif()

find_package(Umfpack)
if(UMFPACK_FOUND)
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${UMFPACK_INCLUDES})
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${UMFPACK_LIBRARIES})
    set(WITH_UMFPACK ON)
else()
    message(WARNING "[Warning] Umfpack not found")
    set(WITH_UMFPACK OFF)
endif()


# set-up module
if(BLAS_FOUND OR OPEN_BLAS_FOUND)
    set(UTOPIA_BLAS TRUE PARENT_SCOPE)
    set(WITH_BLAS ON PARENT_SCOPE)

    list(APPEND BLAS_MODULES
        .
        solvers
        eigensolvers
        eval
    )

    if(LAPACK_FOUND)
        set(WITH_LAPACK ON PARENT_SCOPE)
        set(UTOPIA_LAPACK TRUE PARENT_SCOPE)
        list(APPEND BLAS_MODULES solvers/linear/lapack eigensolvers/lapack)
    endif()

    if(UMFPACK_FOUND)
        set(WITH_UMFPACK ON PARENT_SCOPE)
        set(UTOPIA_UMFPACK TRUE PARENT_SCOPE)
        list(APPEND BLAS_MODULES solvers/linear/umfpack)
    endif()


    utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${BLAS_MODULES}")
    target_include_directories(utopia PUBLIC ${BLAS_INCLUDE_DIR})
    target_link_libraries(utopia PUBLIC ${BLAS_LIBRARIES})
    target_compile_definitions(utopia PUBLIC ${BLAS_DEFINITIONS})

    if(LAPACK_FOUND)
        target_include_directories(utopia PUBLIC ${LAPACK_INCLUDE_DIR})
        target_link_libraries(utopia PUBLIC ${LAPACK_LIBRARIES})
        target_compile_definitions(utopia PUBLIC ${LAPACK_DEFINITIONS})
    endif()

    if(UMFPACK_FOUND)
        target_include_directories(utopia PUBLIC ${UMFPACK_INCLUDES})
        target_link_libraries(utopia PUBLIC ${UMFPACK_LIBRARIES})
    endif()

    set(UTOPIA_THIRDPARTY_INCLUDES ${UTOPIA_THIRDPARTY_INCLUDES} PARENT_SCOPE)
    set(UTOPIA_THIRDPARTY_LIBRARIES ${UTOPIA_THIRDPARTY_LIBRARIES} PARENT_SCOPE)
    set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
endif()
