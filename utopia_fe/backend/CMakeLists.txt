# add all available backends
add_subdirectory(libmesh)
add_subdirectory(homemade)
add_subdirectory(intrepid2)
add_subdirectory(moonolith)
add_subdirectory(stk)
add_subdirectory(arborx)

set(UTOPIA_FE_LIBRARIES
    ${UTOPIA_FE_LIBRARIES}
    PARENT_SCOPE)
