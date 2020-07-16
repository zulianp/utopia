# find dependencies

find_package(OpenCL)

if(OPENCL_FOUND)
    list(APPEND UTOPIA_LIBRARIES ${OPENCL_LIBRARIES})
    list(APPEND UTOPIA_INCLUDES ${OPENCL_INCLUDE_DIRS})
    set(UTOPIA_OPENCL TRUE PARENT_SCOPE)
    set(WITH_OPENCL ON PARENT_SCOPE)
    set(WITH_UTOPIA_OPENCL ON PARENT_SCOPE)
else()
    set(OPENCL_FOUND FALSE)
endif()

# set-up module
if(OPENCL_FOUND)
    utopia_add_library(utopia_opencl ${CMAKE_CURRENT_SOURCE_DIR} "." "utopia_core;utopia_utils;utopia_solvers")
    target_include_directories(utopia_opencl PUBLIC ${OPENCL_INCLUDE_DIRS})
    target_link_libraries(utopia_opencl PUBLIC ${OPENCL_LIBRARIES})

    set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
    set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
endif()
