# find dependencies

find_package(OpenCL)

if(OPENCL_FOUND)
    list(APPEND UTOPIA_BUILD_INCLUDES ${OPENCL_LIBRARIES})
    list(APPEND UTOPIA_BUILD_INCLUDES ${OPENCL_INCLUDE_DIRS})
    set(UTOPIA_OPENCL
        TRUE
        PARENT_SCOPE)
    set(WITH_OPENCL
        ON
        PARENT_SCOPE)
    set(UTOPIA_WITH_OPENCL
        ON
        PARENT_SCOPE)
else()
    set(OPENCL_FOUND FALSE)
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

# set-up module
# if(OPENCL_FOUND)
#     utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
#     target_include_directories(utopia PUBLIC ${OPENCL_INCLUDE_DIRS})
#     target_link_libraries(utopia PUBLIC ${OPENCL_LIBRARIES})

#     set(UTOPIA_BUILD_INCLUDES
#         ${UTOPIA_BUILD_INCLUDES}
#         PARENT_SCOPE)
#     set(UTOPIA_DEP_LIBRARIES
#         ${UTOPIA_DEP_LIBRARIES}
#         PARENT_SCOPE)
# endif()
