list(APPEND SCRIPTING_MODULES
    .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR SCRIPTING_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_c SHARED ${LOCAL_SOURCES})
target_link_libraries(utopia_c utopia)
target_include_directories(utopia_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS utopia_c
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


################ SWIG ##################

find_package(SWIG)

if(SWIG_FOUND)
    include(${SWIG_USE_FILE})

    find_package(PythonLibs)

    if(PythonLibs_FOUND)
        include_directories(${PYTHON_INCLUDE_PATH})
        include_directories(${CMAKE_CURRENT_SOURCE_DIR})
        # set(CMAKE_SWIG_FLAGS "")

        set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
        # set_property(SOURCE utopia.i PROPERTY USE_LIBRARY_INCLUDE_DIRECTORIES TRUE)
        set_source_files_properties(utopia.i PROPERTIES SWIG_FLAGS "-includeall")

        swig_add_library(utopia TYPE SHARED LANGUAGE python SOURCES utopia.i utopia_script.cpp)

        # set_property(TARGET _utopia PROPERTY SWIG_USE_LIBRARY_INCLUDE_DIRECTORIES TRUE)

        swig_link_libraries(utopia ${PYTHON_LIBRARIES} utopia)
        set_target_properties(_utopia PROPERTIES SUFFIX ".so")

        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

        install(
          TARGETS ${SWIG_MODULE_utopia_REAL_NAME}
          # ... add other arguments to install() as necessary
          LIBRARY DESTINATION lib
        )

    endif()

endif()
