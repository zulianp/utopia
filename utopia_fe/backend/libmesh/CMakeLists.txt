utopiafe_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
find_package(LIBMESH)

if(LIBMESH_FOUND)
    target_include_directories(utopia_fe BEFORE PUBLIC ${LIBMESH_INCLUDE_DIRS} ${LIBMESH_INCLUDE_DIR})
    target_link_libraries(utopia_fe INTERFACE ${LIBMESH_LIBRARIES})
    target_include_directories(utopia_fe PRIVATE ${PROJECT_SOURCE_DIR}/rapidxml)

    set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIRS} PARENT_SCOPE)
    set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARIES} PARENT_SCOPE)
endif()
