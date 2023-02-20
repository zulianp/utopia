# find dependencies ################################################### 
if(PITZ_DORA)
    find_package(OpenBLAS)
    if(OPEN_BLAS_FOUND)
        list(APPEND UTOPIA_BUILD_INCLUDES ${BLAS_INCLUDE_DIR})
        list(APPEND UTOPIA_DEP_LIBRARIES ${BLAS_LIBRARIES})
        list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
        set(UTOPIA_BLAS TRUE PARENT_SCOPE)
        set(UTOPIA_WITH_OPEN_BLAS ON)
        # set(UTOPIA_BLAS_DIR ${BLAS_blas_LIBRARY} PARENT_SCOPE)
        
        # # BLAS_blas_LIBRARY
        # set(UTOPIA_BLAS_DIR ${BLAS_VERSION} PARENT_SCOPE)
    else()
        message(WARNING "[Warning] blas not found")
        set(UTOPIA_WITH_BLAS OFF)
    endif()
else()
    SET(OPEN_BLAS_FOUND FALSE)
    find_package(BLAS)
    if(BLAS_FOUND)
        list(APPEND UTOPIA_BUILD_INCLUDES ${BLAS_INCLUDE_DIR})
        list(APPEND UTOPIA_DEP_LIBRARIES ${BLAS_LIBRARIES})
        list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
        set(UTOPIA_WITH_BLAS ON)
        set(UTOPIA_BLAS TRUE PARENT_SCOPE)
    else()
        SET(BLAS_FOUND FALSE)
    endif()
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${LAPACK_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${LAPACK_DEFINITIONS})
    set(UTOPIA_WITH_LAPACK ON)
    set(UTOPIA_LAPACK TRUE PARENT_SCOPE)
else()
    message(WARNING "[Warning] lapack not found")
    set(UTOPIA_WITH_LAPACK OFF)
endif()

find_package(UMFPACK)
if(UMFPACK_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${UMFPACK_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${UMFPACK_LIBRARIES})
    set(UTOPIA_WITH_UMFPACK ON)
    set(UTOPIA_UMFPACK TRUE PARENT_SCOPE)           
else()
    message(WARNING "[Warning] Umfpack not found")
    set(UTOPIA_WITH_UMFPACK OFF)
endif()


# find dependencies ###################################################

if(BLAS_FOUND OR OPEN_BLAS_FOUND)
list(APPEND BLAS_MODULES
        .
        solvers
        eigensolvers
        eval
    )

if(LAPACK_FOUND)
    list(APPEND BLAS_MODULES solvers/linear/lapack eigensolvers/lapack)
endif()

if(UMFPACK_FOUND)
    list(APPEND BLAS_MODULES solvers/linear/umfpack)
endif()


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${BLAS_MODULES}" UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)


set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)


endif()




# # set-up module
# if(BLAS_FOUND OR OPEN_BLAS_FOUND)
#     set(UTOPIA_BLAS TRUE PARENT_SCOPE)
#     set(UTOPIA_WITH_BLAS ON PARENT_SCOPE)

#     list(APPEND BLAS_MODULES
#         .
#         solvers
#         eigensolvers
#         eval
#     )

#     if(LAPACK_FOUND)
#         set(UTOPIA_WITH_LAPACK ON PARENT_SCOPE)
#         set(UTOPIA_LAPACK TRUE PARENT_SCOPE)
#         list(APPEND BLAS_MODULES solvers/linear/lapack eigensolvers/lapack)
#     endif()

#     if(UMFPACK_FOUND)
#         set(UTOPIA_WITH_UMFPACK ON PARENT_SCOPE)
#         set(UTOPIA_UMFPACK TRUE PARENT_SCOPE)
#         list(APPEND BLAS_MODULES solvers/linear/umfpack)
#     endif()


#     utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${BLAS_MODULES}")
#     target_include_directories(utopia PUBLIC ${BLAS_INCLUDE_DIR})
#     target_link_libraries(utopia PUBLIC ${BLAS_LIBRARIES})
#     target_compile_definitions(utopia PUBLIC ${BLAS_DEFINITIONS})

#     if(LAPACK_FOUND)
#         target_include_directories(utopia PUBLIC ${LAPACK_INCLUDE_DIR})
#         target_link_libraries(utopia PUBLIC ${LAPACK_LIBRARIES})
#         target_compile_definitions(utopia PUBLIC ${LAPACK_DEFINITIONS})
#     endif()

#     if(UMFPACK_FOUND)
#         target_include_directories(utopia PUBLIC ${UMFPACK_INCLUDES})
#         target_link_libraries(utopia PUBLIC ${UMFPACK_LIBRARIES})
#     endif()

#     set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES} PARENT_SCOPE)
#     set(UTOPIA_DEP_LIBRARIES ${UTOPIA_DEP_LIBRARIES} PARENT_SCOPE)
#     set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
# endif()







