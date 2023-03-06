find_package(yaml-cpp HINTS ${CMAKE_SOURCE_DIR}/../../dependencies/ ${CMAKE_SOURCE_DIR}/../../dependencies/yaml-cpp/share/cmake/yaml-cpp)

if(yaml-cpp_FOUND)
    set(UTOPIA_ENABLE_YAML_CPP ON PARENT_SCOPE)
    set(UTOPIA_ENABLE_YAML_CPP TRUE PARENT_SCOPE)

    set(UTOPIA_YAML_CPP_DIR ${yaml-cpp_DIR} PARENT_SCOPE)
    

    list(
        APPEND
        YAMLCPP_MODULES
        .
        )

    # utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${YAMLCPP_MODULES}")

    get_target_property(YAML_CPP_INCLUDE_DIR yaml-cpp
                        INTERFACE_INCLUDE_DIRECTORIES)

    get_target_property(YAML_CPP_LIBRARIES yaml-cpp
                        IMPORTED_LOCATION)

    if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp
                            IMPORTED_LOCATION_RELEASE)
    endif()
    if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp
                            IMPORTED_LOCATION_DEBUG)
    endif()
    if(NOT YAML_CPP_LIBRARIES)
        get_target_property(YAML_CPP_LIBRARIES yaml-cpp
                            IMPORTED_LOCATION_NOCONFIG)
    endif()

    message(STATUS "yaml-cpp found! Includes: ${YAML_CPP_INCLUDE_DIR}\nLibrary: ${YAML_CPP_LIBRARIES}\n")

    # target_link_libraries(utopia PUBLIC yaml-cpp)
    # target_include_directories(utopia PUBLIC ${YAML_CPP_INCLUDE_DIR})
    # target_link_libraries(utopia PUBLIC ${YAML_CPP_LIBRARIES})

    list(APPEND UTOPIA_BUILD_INCLUDES ${YAML_CPP_INCLUDE_DIR})
    list(APPEND UTOPIA_DEP_LIBRARIES ${YAML_CPP_LIBRARIES})

    set(UTOPIA_BUILD_INCLUDES
        ${UTOPIA_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_DEP_LIBRARIES
        ${UTOPIA_DEP_LIBRARIES}
        PARENT_SCOPE)

    scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${YAMLCPP_MODULES}" UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)


set(UTOPIA_BUILD_INCLUDES
    ${UTOPIA_BUILD_INCLUDES}
    PARENT_SCOPE)

set(UTOPIA_HEADERS
    ${UTOPIA_HEADERS}
    PARENT_SCOPE)

set(UTOPIA_SOURCES
    ${UTOPIA_SOURCES}
    PARENT_SCOPE)

else()
    include(../../cmake/InstallYAMLCPP.cmake)
endif()