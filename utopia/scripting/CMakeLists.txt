# utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} ".")


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "."
                 UTOPIA_DEP_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

set(UTOPIA_DEP_INCLUDES
    ${UTOPIA_DEP_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)


# ############### SWIG ##################
cmake_policy(SET "CMP0086" NEW)

find_package(SWIG REQUIRED)

if(SWIG_FOUND)
    include(${SWIG_USE_FILE})

    find_package(Python3 COMPONENTS Development NumPy REQUIRED)

    if(Python3_FOUND)
        set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
        set_source_files_properties(utopia.i PROPERTIES SWIG_FLAGS
                                                        "-includeall")

        set_source_files_properties(SOURCE utopia.i PROPERTY
                                    SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)

        message(STATUS "Found python libraries: ${Python_LIBRARIES}")
        swig_add_library(
            utopya
            TYPE SHARED
            LANGUAGE python
            SOURCES utopia.i utopia_script.hpp utopia_script.cpp)

        # swig_link_libraries(utopya PUBLIC Python3::Python Python3::NumPy utopia)
        swig_link_libraries(utopya PUBLIC Python3::Module utopia)
        set_target_properties(utopya PROPERTIES SUFFIX ".so")

        set_property(TARGET utopya PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES
                                            TRUE)

        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia.py DESTINATION lib)

        install(TARGETS ${SWIG_MODULE_utopya_REAL_NAME} LIBRARY DESTINATION lib)

        # find_package(NumPy REQUIRED)
        message(STATUS "PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")

    else()
        message(WARNING "Could not find Python - Disabling scripting support")
    endif()
else()
    message(WARNING "Could not find SWIG - Disabling scripting support")
endif()
