find_package(M3ELINSOL)

if(M3ELINSOL_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${M3ELINSOL_INCLUDES})
    list(APPEND UTOPIA_DEP_LIBRARIES ${M3ELINSOL_LIBRARIES})
    list(APPEND UTOPIA_DEFS ${M3ELINSOL_DEFS})
    set(UTOPIA_WITH_M3ELINSOL ON)
else()
    message(WARNING "[Warning] M3ELINSOL not found")
    set(UTOPIA_WITH_M3ELINSOL OFF)
endif()



scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "." UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)


set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)



# # set-up module
# if(M3ELINSOL_FOUND)
#     set(UTOPIA_M3ELINSOL TRUE PARENT_SCOPE)
#     set(UTOPIA_WITH_M3ELINSOL ON PARENT_SCOPE)


#     utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
#     target_sources(utopia PRIVATE ${M3ELINSOL_SOURCE_FILES})
#     target_include_directories(utopia PUBLIC ${M3ELINSOL_INCLUDES})
#     target_link_libraries(utopia PUBLIC ${M3ELINSOL_LIBRARIES})

#     find_package(OpenMP)
#     if(OPENMP_FOUND)
#         SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" PARENT_SCOPE)
#         SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" PARENT_SCOPE)
#         SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" PARENT_SCOPE)
#         set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS} -lgfortran" PARENT_SCOPE)
#     endif()

#     # find_package(gfortran)
#     # if(GFORTRAN_FOUND)
#     # 	message(STATUS "Found gfortran")
#     # endif()

#     # set_target_properties(utopia_m3elinsol PROPERTIES COMPILE_FLAGS "-fopenmp")
#     # set_target_properties(utopia_m3elinsol PROPERTIES LINK_FLAGS "-fopenmp")

#     if(TARGET utopia_petsc)
#         target_link_libraries(utopia PUBLIC utopia_petsc)
#     endif()

#     if(TARGET utopia_blas)
#         target_link_libraries(utopia PUBLIC utopia_blas)
#     endif()

#     if(TARGET utopia_trilinos)
#         target_link_libraries(utopia PUBLIC utopia_trilinos)
#     endif()

#     target_compile_definitions(utopia PUBLIC ${M3ELINSOL_DEFS})

#     set(UTOPIA_BUILD_INCLUDES ${UTOPIA_BUILD_INCLUDES} PARENT_SCOPE)
#     set(UTOPIA_DEP_LIBRARIES ${UTOPIA_DEP_LIBRARIES} PARENT_SCOPE)
#     set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
# endif()
