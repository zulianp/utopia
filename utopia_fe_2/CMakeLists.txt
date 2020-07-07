cmake_minimum_required(VERSION 3.8.1)
project(utopia_fe)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

#To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)

option(ENABLE_SANITIZER "check for memory access problems" OFF)
option(ENABLE_MOONOLITH_PROFILER "enable moonolith profiling capabilities" OFF)
option(ENABLE_CXX14_FEATURES "Enable usage of cxx14 standard" ON)
option(ENABLE_MARS "Tries to use mars backend" ON)
option(ENABLE_INTREPID2 "Tries to use the intrepid2 related features" OFF)
option(ENABLE_LIBMESH "Tries to use the LibMesh related features" ON)
option(ENABLE_MOONOLITH "enable moonolith transfer capabilities" ON)
#################################################################################
#################################################################################
#################################################################################

set(UTOPIA_FE_DEV_FLAGS
        "-Wall -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder"
)

set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_FE_DEV_FLAGS} -Werror=unused-variable -Werror=unused-local-typedef -Werror=deprecated-declarations ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")


if(UNIX AND NOT APPLE)
       set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(ENABLE_TRILINOS_ALGEBRA)
    set(WITH_TRILINOS_ALGEBRA TRUE)
endif()

if(ENABLE_NEW_TRANSFER)
    set(WITH_NEW_TRANSFER TRUE)
endif()

set(UTOPIA_FE_ROOT_PATH ${CMAKE_SOURCE_DIR})
set(UTOPIA_FE_CMAKES_PATH ${UTOPIA_FE_ROOT_PATH}/libcmake)


if (ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
        SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ")
    endif()
endif()

################################## Include cmake functions and find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_FE_ROOT_PATH}/cmake")
include(${UTOPIA_FE_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

file(GLOB_RECURSE FUNC_LIB "{UTOPIA_FE_ROOT_PATH}/cmake/utils" )
foreach(funcLib ${FUNC_LIB})
    include(${funcLib})
endforeach(funcLib)

find_package(Utopia REQUIRED)
add_definitions(${UTOPIA_DEFS})
include_directories(${UTOPIA_INCLUDES})
include_directories(rapidxml)


########################################################################
function(find_project_files rootPath dirPaths headers sources)
    SET(verbose TRUE)

    SET(theaders ${${headers}})
    SET(tsources ${${sources}})

    SET(ROOT_PATH ${${rootPath}})
    SET(DIR_PATHS ${${dirPaths}})

    FOREACH(INCLUDE_PATH ${DIR_PATHS})
        INCLUDE_DIRECTORIES(BEFORE ${ROOT_PATH}/${INCLUDE_PATH})

        file(GLOB TEMP_HPPSRC CONFIGURE_DEPENDS "${ROOT_PATH}/${INCLUDE_PATH}/*.cpp" )
        file(GLOB TEMP_SRC CONFIGURE_DEPENDS "${ROOT_PATH}/${INCLUDE_PATH}/*.c" )
        file(GLOB TEMP_HPPDR CONFIGURE_DEPENDS "${ROOT_PATH}/${INCLUDE_PATH}/*.hpp" )
        file(GLOB TEMP_HDR CONFIGURE_DEPENDS "${ROOT_PATH}/${INCLUDE_PATH}/*.h" )

        SOURCE_GROUP(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})

        SET(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
        SET(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
    ENDFOREACH(INCLUDE_PATH)

    SET(${headers} ${theaders} PARENT_SCOPE)
    SET(${sources} ${tsources} PARENT_SCOPE)
endfunction()
########################################################################

# all arguments after module_name are the dependencies, we access them via ${ARGN}
macro(add_utopiafe_module module_name)
    add_library(${module_name} STATIC ${LOCAL_SOURCES})
    target_link_libraries(${module_name} PUBLIC ${ARGN} ${UTOPIA_LIBRARIES})
    target_include_directories(${module_name} BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
    foreach(module_dependency ${ARGN})
        target_include_directories(${module_name} PUBLIC $<TARGET_PROPERTY:${module_dependency},INTERFACE_INCLUDE_DIRECTORIES>)
    endforeach()
    target_include_directories(${module_name} PRIVATE ${UTOPIA_INCLUDES})
    set_utopia_compiler_features(${module_name})
        if (NOT TOP_LEVEL)
        set(UTOPIA_FE_LIBRARIES "${UTOPIA_FE_LIBRARIES};${module_name}" PARENT_SCOPE)
    endif()
    install(TARGETS ${module_name}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
    install(FILES ${LOCAL_HEADERS} DESTINATION include)
endmacro(add_utopiafe_module)


macro(add_utopiafe_app app_name)
    add_library(${app_name} EXCLUDE_FROM_ALL STATIC ${LOCAL_SOURCES})
    target_link_libraries(${app_name} PUBLIC ${ARGN} ${UTOPIA_LIBRARIES})
    target_include_directories(${app_name} BEFORE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
    foreach(app_dependency ${ARGN})
        target_include_directories(${app_name} PUBLIC $<TARGET_PROPERTY:${app_dependency},INTERFACE_INCLUDE_DIRECTORIES>)
    endforeach()
    target_include_directories(${app_name} PRIVATE ${UTOPIA_INCLUDES})
    set_utopia_compiler_features(${app_name})
        if (NOT TOP_LEVEL)
        set(UTOPIA_FE_APPS "${UTOPIA_FE_APPS};${app_name}" PARENT_SCOPE)
    endif()
    # install(TARGETS ${app_name}
    #     RUNTIME DESTINATION bin
    #     LIBRARY DESTINATION lib
    #     ARCHIVE DESTINATION lib
    # )
    # install(FILES ${LOCAL_HEADERS} DESTINATION include)
endmacro(add_utopiafe_app)

list(APPEND UTOPIA_FE_MODULES
    backend/libmesh
    backend/intrepid2
    backend/mars
    backend/moonolith
    backend/petsc
    kernels
    core
    base
 )

if(INSTALL_LIBMESH)
    include(InstallLIBMESH)
endif()

if(ENABLE_INTREPID2)
    list(APPEND UTOPIA_FE_MODULES backend/intrepid2)
endif()

include_directories(BEFORE ${CMAKE_SOURCE_DIR})
include_directories(BEFORE ${CMAKE_BINARY_DIR})

FOREACH(MODULE ${UTOPIA_FE_MODULES})
        add_subdirectory("${MODULE}")
ENDFOREACH(MODULE)

configure_file(${UTOPIA_FE_ROOT_PATH}/utopia_fe_config.hpp.in ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.hpp DESTINATION include)

#create main library
list(APPEND LOCAL_HEADERS utopia_fe.hpp utopia_fe_Instance.hpp)
list(APPEND LOCAL_SOURCES utopia_fe.cpp utopia_fe_Instance.cpp)
set(TOP_LEVEL 1)
add_utopiafe_module(utopia_fe ${UTOPIA_FE_LIBRARIES})
add_executable(utopia_fe_exec EXCLUDE_FROM_ALL main.cpp ${LOCAL_HEADERS})
target_include_directories(utopia_fe_exec PRIVATE ${UTOPIA_INCLUDES})
target_link_libraries(utopia_fe_exec utopia_fe ${UTOPIA_FE_APPS})

include(apps/apps.cmake)

install(TARGETS utopia_fe
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

function(append_list_to_string_for_makefile OUTPUT_STRING PRE_LIST POST_LIST ITEM_PREFIX LIST_OF_STRINGS)
    #http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

    LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

    SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    FOREACH(LIST_EL ${LIST_OF_STRINGS})

        string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
        if(IS_FRAME_WORK)
            get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
            string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
            string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")

            SET(TEMPSTR "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
        else()
            SET(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL} ")
        endif()
    ENDFOREACH(LIST_EL)

    SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()

function(append_list_to_string_for_makefile_with_postfix OUTPUT_STRING PRE_LIST POST_LIST ITEM_PREFIX ITEM_POSTFIX LIST_OF_STRINGS)
    #http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

    LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

    SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    FOREACH(LIST_EL ${LIST_OF_STRINGS})

        string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
        if(IS_FRAME_WORK)
            get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
            string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
            string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")

            SET(TEMPSTR "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
        else()
            SET(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL}${ITEM_POSTFIX} ")
        endif()
    ENDFOREACH(LIST_EL)

    SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()

#################################################################################

function(append_list_to_string_for_cmake OUTPUT_STRING PRE_LIST POST_LIST LIST_OF_STRINGS)
#http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro
    LIST(REMOVE_DUPLICATES LIST_OF_STRINGS)

    SET(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    FOREACH(LIST_EL ${LIST_OF_STRINGS})
        SET(TEMPSTR "${TEMPSTR}\t\"${LIST_EL}\"\n")
    ENDFOREACH(LIST_EL)

    SET(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    SET(${OUTPUT_STRING} ${TEMPSTR} PARENT_SCOPE)
endfunction()



# Generate files useful for dependent build systems
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n#Includes the content of utopia-config.makefile generated by utopia\n\n")

list(APPEND UTOPIA_FE_INCLUDES ${UTOPIA_INCLUDES})
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include -I${LIBMESH_INCLUDE_DIR} " " " "-I" "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" "\n" "-I" "${MOONOLITH_INCLUDES}")

append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib -lutopia_fe " " " "-l" "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " "" "${LIBMESH_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " "" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "" " " "" "${MOONOLITH_LIBRARIES}")

set(CONFIG_CONTENT "${CONFIG_CONTENT}\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_FE_DEFS=" "\n" "" "${UTOPIA_DEFS}")


set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_FLAGS=${UTOPIA_CXX_FLAGS}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_CXX_COMPILER=${UTOPIA_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_FE_C_COMPILER=${UTOPIA_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile ${CONFIG_CONTENT})



set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n"        "" "${UTOPIA_FE_INCLUDES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "" "" "${LIBMESH_INCLUDE_DIR}")

if(MOONOLITH_INCLUDES)
append_list_to_string_for_cmake(CONFIG_CONTENT "" ")\n\n"  ${MOONOLITH_INCLUDES})
endif()

append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib" "" "\n\t\"-l" "\"" "${UTOPIA_FE_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "\n\t-L${CMAKE_INSTALL_PREFIX}/lib" "" "\n\t\"" "\"" "${MOONOLITH_LIBRARIES}")
append_list_to_string_for_makefile_with_postfix(CONFIG_CONTENT "" "" "\n\t\"" "\"" "${LIBMESH_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" "\n)\n\n" "${UTOPIA_LIBRARIES}")

append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_DEFS\n"           "" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT ""           "\n)\n\n" "${UTOPIA_FE_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_FLAGS\n"      "\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_C_COMPILER\n"     "\n)\n\n" "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_FE_CXX_COMPILER\n"   "\n)\n\n" "${CMAKE_CXX_COMPILER}")


file(WRITE ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake ${CONFIG_CONTENT})

install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia_fe_config.cmake    DESTINATION config)

##########################################################################################

# test suite
add_executable(utopia_fe_test EXCLUDE_FROM_ALL test.cpp)
target_link_libraries(utopia_fe_test utopia_fe)
include(tests/tests.cmake)

##########################################################################################

add_custom_target(apps)
add_dependencies(apps utopia_fe utopia_fe_exec)

##########################################################################################

add_custom_target(complete)
add_dependencies(complete apps utopia_fe_test)

##########################################################################################

