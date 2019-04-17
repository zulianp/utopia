# find dependencies

find_package(OpenCL)

if(OPENCL_FOUND)
    if(OPENCL_FOUND)
        set(OPENCL_FOUND TRUE)
        list(APPEND UTOPIA_LIBRARIES ${OPENCL_LIBRARIES})
        list(APPEND UTOPIA_INCLUDES ${OPENCL_INCLUDE_DIRS})
    else()
        set(OPENCL_FOUND FALSE)
    endif()

endif()

# set-up module
if(OPENCL_FOUND)
    set(UTOPIA_OPENCL TRUE PARENT_SCOPE)
    set(WITH_OPENCL ON PARENT_SCOPE)
    set(WITH_UTOPIA_OPENCL ON PARENT_SCOPE)

    list(APPEND OPENCL_MODULES
        .
        )

    find_project_files(CMAKE_CURRENT_SOURCE_DIR OPENCL_MODULES LOCAL_HEADERS LOCAL_SOURCES)

    add_library(utopia_opencl STATIC utopia_opencl.cpp ${LOCAL_SOURCES})

    foreach(MODULE ${OPENCL_MODULES})
        target_include_directories(utopia_opencl PUBLIC ${MODULE})
    endforeach(MODULE)

    target_include_directories(utopia_opencl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_include_directories(utopia_opencl PUBLIC ${OPENCL_INCLUDE_DIRS})

    target_link_libraries(utopia_opencl utopia_core utopia_utils utopia_solvers)
    target_link_libraries(utopia_opencl ${OPENCL_LIBRARIES})


    install(TARGETS utopia_opencl
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

    install(FILES ${LOCAL_HEADERS} DESTINATION include)

    set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_opencl" PARENT_SCOPE)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")


    set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
    set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
    set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
endif()
