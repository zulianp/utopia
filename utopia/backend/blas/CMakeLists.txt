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
#     set(UTOPIA_ENABLE_BLAS TRUE PARENT_SCOPE)
#     set(UTOPIA_ENABLE_BLAS ON PARENT_SCOPE)

#     list(APPEND BLAS_MODULES
#         .
#         solvers
#         eigensolvers
#         eval
#     )

#     if(LAPACK_FOUND)
#         set(UTOPIA_ENABLE_LAPACK ON PARENT_SCOPE)
#         set(UTOPIA_ENABLE_LAPACK TRUE PARENT_SCOPE)
#         list(APPEND BLAS_MODULES solvers/linear/lapack eigensolvers/lapack)
#     endif()

#     if(UMFPACK_FOUND)
#         set(UTOPIA_ENABLE_UMFPACK ON PARENT_SCOPE)
#         set(UTOPIA_ENABLE_UMFPACK TRUE PARENT_SCOPE)
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







