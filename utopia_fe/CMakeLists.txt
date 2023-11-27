cmake_minimum_required(VERSION 3.13)
# project(utopia_fe)
cmake_policy(SET CMP0079 NEW)
project(
    utopia_fe
    VERSION 0.1.0
    LANGUAGES "CXX" "C"
    HOMEPAGE_URL "https://bitbucket.org/zulianp/utopia"
    DESCRIPTION
        "Library for portable finite element discretizations."
)


# ##############################################################################
# ######################## UTOPIA_FE DEFAULT INCLUDES
# ##############################################################################
# ##############################################################################

include(cmake/UtopiaFEOptions.cmake)
include(cmake/UtopiaFEFunctions.cmake)
include(cmake/UtopiaFEDefaults.cmake)


# ##############################################################################
# ##############################################################################
# ##############################################################################

if(NOT UTOPIA_ENABLE_MARS)
    if(UTOPIA_ENABLE_WARNINGS)

        set(UTOPIA_FE_DEV_FLAGS
            "-Wall -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder -Werror=return-type" # -Werror=uninitialized
)

        set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Wextra ")

    else()
        set(UTOPIA_FE_DEV_FLAGS "-w -Wno-deprecated-declarations")
    endif()

endif()

# More annoying set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_FE_DEV_FLAGS} -Wextra
# -Werror=uninitialized ")

# More restrictive set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_FE_DEV_FLAGS}
# -Werror=unused-variable -Werror=unused-local-typedef
# -Werror=deprecated-declarations ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(UTOPIA_ENABLE_TRILINOS_ALGEBRA)
    set(UTOPIA_ENABLE_TRILINOS_ALGEBRA TRUE)
endif()

if(UTOPIA_ENABLE_NEW_TRANSFER)
    set(UTOPIA_ENABLE_NEW_TRANSFER TRUE)
endif()

set(UTOPIA_FE_ROOT_PATH ${CMAKE_SOURCE_DIR})
set(UTOPIA_FE_CMAKES_PATH ${UTOPIA_FE_ROOT_PATH}/libcmake)

if(UTOPIA_ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
        set(CMAKE_CXX_FLAGS_DEBUG
            "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls "
        )
    endif()
endif()



# ################################# Include cmake functions and find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_FE_ROOT_PATH}/cmake")
include(${UTOPIA_FE_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

if(UTOPIA_ENABLE_FLUYA_MODE)
include(EnableFluyaMode)
endif()


file(GLOB_RECURSE FUNC_LIB "{UTOPIA_FE_ROOT_PATH}/cmake/utils")
foreach(funcLib ${FUNC_LIB})
    include(${funcLib})
endforeach(funcLib)

find_package(Utopia REQUIRED)
add_definitions(${UTOPIA_DEFS})

######################################################################

list(APPEND UTOPIA_FE_MODULES frontend interoperability)

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND UTOPIA_FE_MODULES backend/libmesh edsl)
endif()

if(UTOPIA_ENABLE_MOONOLITH)
    list(APPEND UTOPIA_FE_MODULES backend/moonolith)
endif()

if(UTOPIA_ENABLE_STK)
    list(APPEND UTOPIA_FE_MODULES backend/stk)
endif()

if(UTOPIA_ENABLE_PETSC)
    list(APPEND UTOPIA_FE_MODULES backend/petsc)
endif()


##################################################

if(UTOPIA_ENABLE_INTREPID2 OR UTOPIA_ENABLE_MARS OR UTOPIA_ENABLE_LIBMESH_KOKKOS)
    list(APPEND UTOPIA_FE_MODULES backend/kokkos)
endif()

if(UTOPIA_ENABLE_INTREPID2)
    list(APPEND UTOPIA_FE_MODULES backend/intrepid2)
endif()

if(UTOPIA_ENABLE_MARS)
    list(APPEND UTOPIA_FE_MODULES backend/mars)
endif()

if(UTOPIA_ENABLE_ARBORX)
    list(APPEND UTOPIA_FE_MODULES backend/arborx)
endif()

##################################################


list(APPEND UTOPIA_FE_MODULES base core)

if(UTOPIA_ENABLE_LIBMESH)
    include(InstallLIBMESH)
endif()


add_library(utopia_fe utopia_fe.cpp)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# if(UTOPIA_ENABLE_MARS_ADIOS2)
#     find_package(adios2 REQUIRED)
#     target_link_libraries(utopia_fe PUBLIC adios2::adios2)
# endif()

set_utopia_compiler_features(utopia_fe)
target_include_directories(
    utopia_fe BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                            $<INSTALL_INTERFACE:include>)
target_include_directories(
    utopia_fe BEFORE PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(
    TARGETS utopia_fe
    EXPORT UtopiaFETargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
target_link_libraries(utopia_fe PUBLIC Utopia::utopia)
install(FILES utopia_fe.hpp DESTINATION include)

# now add all subdirectories to utopia_fe library
foreach(MODULE ${UTOPIA_FE_MODULES})
    add_subdirectory("${MODULE}")
endforeach(MODULE)

configure_file(${UTOPIA_FE_ROOT_PATH}/utopia_fe_config.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp DESTINATION include)

set(TOP_LEVEL 1)
add_executable(utopia_fe_exec EXCLUDE_FROM_ALL main.cpp ${LOCAL_HEADERS})
# target_link_libraries(utopia_fe_exec PRIVATE utopia_fe ${UTOPIA_FE_APPS})

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------


if(UTOPIA_FE_LIBRARIES_TRAILING)
    list(APPEND UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES_TRAILING})
endif()

# Generate files useful for dependent build systems
set(CONFIG_CONTENT
    "#Automatically generated by utopia CMakeLists.txt\n#Includes the content of utopia-config.makefile generated by utopia\n\n"
)

list(APPEND UTOPIA_FE_INCLUDES ${UTOPIA_INCLUDES})
append_list_to_string_for_makefile(
    CONFIG_CONTENT
    "UTOPIA_FE_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include -I${LIBMESH_INCLUDE_DIR} "
    " "
    "-I"
    "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" "\n" "-I"
                                   "${MOONOLITH_INCLUDES}")

message("MOONOLITH_INCLUDES=${MOONOLITH_INCLUDES}")

append_list_to_string_for_makefile(
    CONFIG_CONTENT
    "UTOPIA_FE_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib -lutopia_fe " " " "-l"
    "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${LIBMESH_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " ""
                                   "${MOONOLITH_LIBRARIES}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_DEFS=" "\n" ""
                                   "${UTOPIA_DEFS}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_FLAGS=${UTOPIA_CXX_FLAGS}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_FE_CXX_COMPILER=${UTOPIA_CXX_COMPILER}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_FE_C_COMPILER=${UTOPIA_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(
    CONFIG_CONTENT
    "LIST(APPEND UTOPIA_FE_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n" ""
    "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "" "${LIBMESH_INCLUDE_DIR}")

if(MOONOLITH_INCLUDES)
    append_list_to_string_for_cmake(CONFIG_CONTENT "" ")\n\n"
                                    ${MOONOLITH_INCLUDES})
else()
    set(CONFIG_CONTENT "${CONFIG_CONTENT})")
endif()

append_list_to_string_for_makefile_with_postfix(
    CONFIG_CONTENT
    "LIST(APPEND UTOPIA_FE_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib\n\t-lutopia_fe"
    ""
    "\n\t\"-l"
    "\""
    "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\""
                                                "\"" "${MOONOLITH_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\""
                                                "\"" "${LIBMESH_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" "\n)\n\n"
                                "${UTOPIA_LIBRARIES}")

append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_DEFS\n"
                                "" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "\n)\n\n" "${UTOPIA_FE_DEFS}")
append_list_to_string_for_cmake(
    CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_FLAGS\n" "\n)\n\n"
    "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(
    CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_C_COMPILER\n" "\n)\n\n"
    "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(
    CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_COMPILER\n" "\n)\n\n"
    "${CMAKE_CXX_COMPILER}")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake ${CONFIG_CONTENT})
list(JOIN CONFIG_CONTENT "\n" OLD_CMAKE_CONFIG_CONTENT)

install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake DESTINATION config)

# new style exported targets, this allows to use target_link_libraries(my_target
# UtopiaFE::utopia_fe)
install(
    EXPORT UtopiaFETargets
    FILE UtopiaFETargets.cmake
    NAMESPACE UtopiaFE::
    DESTINATION lib/cmake)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtopiaFEConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfig.cmake"
    INSTALL_DESTINATION "lib/cmake")
# generate the version file for the config file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfigVersion.cmake"
    VERSION "${UtopiaFE_VERSION_MAJOR}.${UtopiaFE_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion)
# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UtopiaFEConfig.cmake
        DESTINATION lib/cmake/)

add_executable(utopia_fe_test test.cpp)

target_link_libraries(utopia_fe_test utopia_fe)
target_link_libraries(utopia_fe_exec utopia_fe)

include(tests/tests.cmake)
include(apps/apps.cmake)

# if(UTOPIA_ENABLE_LIBMESH AND UTOPIA_ENABLE_LIBMESH_DEPRECATED)
#     include(deprecated/apps/apps.cmake)
# endif()

# ##############################################################################

add_custom_target(apps)
add_dependencies(apps utopia_fe utopia_fe_exec)

# ##############################################################################

add_custom_target(complete)
add_dependencies(complete apps utopia_fe_test)

# ##############################################################################

# include(cmake/PostInstall.cmake)

message(STATUS "UTOPIA_FE_LIBRARIES=${UTOPIA_FE_LIBRARIES}")
