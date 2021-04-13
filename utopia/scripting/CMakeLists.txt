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

    find_package(Python COMPONENTS Development)

    if(Python_FOUND)
        set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
        set_source_files_properties(utopia.i PROPERTIES SWIG_FLAGS
                                                        "-includeall")

        message(STATUS "Found python libraries: ${Python_LIBRARIES}")
        swig_add_library(
            utopya
            TYPE SHARED
            LANGUAGE python
            SOURCES utopia.i utopia_script.cpp)

        swig_link_libraries(utopya Python::Python utopia)
        set_target_properties(utopya PROPERTIES SUFFIX ".so")

        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

        install(TARGETS ${SWIG_MODULE_utopya_REAL_NAME} LIBRARY DESTINATION lib)
    else()
        message(WARNING "Could not find Python - Disabling scripting support")
    endif()
else()
    message(WARNING "Could not find SWIG - Disabling scripting support")
endif()
