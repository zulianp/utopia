utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

################ SWIG ##################

find_package(SWIG)

if(SWIG_FOUND)
    # 2020-07-16 AF: Is this even working???
    #                we are using set_target_properties on the target "_utopia", which is never defined???

    include(${SWIG_USE_FILE})

    find_package(PythonLibs)

    if(PythonLibs_FOUND)
        target_include_directories(utopia PRIVATE ${PYTHON_INCLUDE_PATH})
        # set(CMAKE_SWIG_FLAGS "")

        set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
        # set_property(SOURCE utopia.i PROPERTY USE_LIBRARY_INCLUDE_DIRECTORIES TRUE)
        set_source_files_properties(utopia.i PROPERTIES SWIG_FLAGS "-includeall")

        swig_add_library(utopia TYPE SHARED LANGUAGE python SOURCES utopia.i utopia_script.cpp)

        # set_property(TARGET _utopia PROPERTY SWIG_USE_LIBRARY_INCLUDE_DIRECTORIES TRUE)

        swig_link_libraries(utopia ${PYTHON_LIBRARIES} utopia)
        set_target_properties(utopia PROPERTIES SUFFIX ".so")

        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

        install(
          TARGETS ${SWIG_MODULE_utopia_REAL_NAME}
          # ... add other arguments to install() as necessary
          LIBRARY DESTINATION lib
        )
    else()
        message(WARNING "Could not find Python")
    endif()
else()
        message(WARNING "Could not find Swig")
endif()
