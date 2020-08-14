utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")

# ############### SWIG ##################

find_package(SWIG)

if(SWIG_FOUND)
    # 2020-07-16 AF: Is this even working??? we are using set_target_properties
    # on the target "_utopia", which is never defined??? PZ: This was necessary
    # for dealing with SWIG quirks. For some reason it automatically added the _
    # in front

    # cmake_policy(CMP0078)

    # cmake_policy(SET "CMP0078" NEW)

    include(${SWIG_USE_FILE})

    find_package(PythonLibs)

    if(PythonLibs_FOUND)
        target_include_directories(utopia PRIVATE ${PYTHON_INCLUDE_PATH})
        set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
        set_source_files_properties(utopia.i PROPERTIES SWIG_FLAGS
                                                        "-includeall")

        message(STATUS "${PYTHON_INCLUDE_PATH}")
        swig_add_library(
            utopya
            TYPE SHARED
            LANGUAGE python
            SOURCES utopia.i utopia_script.cpp utopia_swig_Vector.cpp utopia_swig_Matrix.cpp utopia_swig_Communicator.cpp utopia_swig_Layout.cpp)

        swig_link_libraries(utopya ${PYTHON_LIBRARIES} utopia)
        target_include_directories(utopya PRIVATE ${PYTHON_INCLUDE_PATH})
        set_target_properties(utopya PROPERTIES SUFFIX ".so")

        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

        install(TARGETS ${SWIG_MODULE_utopia_REAL_NAME} LIBRARY DESTINATION lib)
    else()
        message(WARNING "Could not find Python")
    endif()
else()
    message(WARNING "Could not find SWIG")
endif()
