# find dependencies

find_package(OpenCL)

if(OPENCL_FOUND)
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${OPENCL_LIBRARIES})
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${OPENCL_INCLUDE_DIRS})
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

# set-up module
if(OPENCL_FOUND)
    utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_include_directories(utopia PUBLIC ${OPENCL_INCLUDE_DIRS})
    target_link_libraries(utopia PUBLIC ${OPENCL_LIBRARIES})

    set(UTOPIA_THIRDPARTY_INCLUDES
        ${UTOPIA_THIRDPARTY_INCLUDES}
        PARENT_SCOPE)
    set(UTOPIA_THIRDPARTY_LIBRARIES
        ${UTOPIA_THIRDPARTY_LIBRARIES}
        PARENT_SCOPE)
endif()
