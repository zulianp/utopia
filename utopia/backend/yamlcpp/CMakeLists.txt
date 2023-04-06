find_package(yaml-cpp HINTS ${CMAKE_SOURCE_DIR}/../../dependencies/ ${CMAKE_SOURCE_DIR}/../../dependencies/yaml-cpp/share/cmake/yaml-cpp)

if(yaml-cpp_FOUND)
    set(UTOPIA_WITH_YAML_CPP TRUE)
    set(UTOPIA_WITH_YAML_CPP TRUE PARENT_SCOPE)

    list(
        APPEND
        YAMLCPP_MODULES
        .
        )

    utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${YAMLCPP_MODULES}")

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

    if(NOT YAML_CPP_LIBRARIES)
        SET(YAML_CPP_LIBRARIES "yaml-cpp")
    endif()

    message(STATUS "yaml-cpp found! Includes: ${YAML_CPP_INCLUDE_DIR}\nLibrary: ${YAML_CPP_LIBRARIES}\n")

    # target_link_libraries(utopia PUBLIC yaml-cpp)
    target_include_directories(utopia PUBLIC ${YAML_CPP_INCLUDE_DIR})
    target_link_libraries(utopia PUBLIC ${YAML_CPP_LIBRARIES})

    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${YAML_CPP_INCLUDE_DIR})
    list(APPEND UTOPIA_THIRDPARTY_LIBRARIES ${YAML_CPP_LIBRARIES})

    set(UTOPIA_THIRDPARTY_INCLUDES
        ${UTOPIA_THIRDPARTY_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_THIRDPARTY_LIBRARIES
        ${UTOPIA_THIRDPARTY_LIBRARIES}
        PARENT_SCOPE)

else()
    include(../../cmake/InstallYAMLCPP.cmake)
endif()
