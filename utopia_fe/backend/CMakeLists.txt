
#add all available backends
add_subdirectory(libmesh)
add_subdirectory(homemade)
set(UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES} PARENT_SCOPE)
