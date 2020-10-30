cmake_minimum_required(VERSION 3.13)

project(utopia)
cmake_policy(SET "CMP0079" NEW)

set(UTOPIA_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/cmake")

# ##############################################################################
# ######################## UTOPIA VERSION INFORMATION
# ##############################################################################
# ##############################################################################

include(UtopiaGitVersion)
include(UtopiaClangFormat)

set(UTOPIA_VERSION_MAJOR 0)
set(UTOPIA_VERSION_MINOR 1)
set(UTOPIA_VERSION_PATCH 0)
set(UTOPIA_VERSION
    "${UTOPIA_VERSION_MAJOR}.${UTOPIA_VERSION_MINOR}.${UTOPIA_VERSION_PATCH}")
# ##############################################################################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

# To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)

set(CMAKE_CXX_EXTENSIONS OFF)
# ##############################################################################
# ##############################################################################
# ##############################################################################

# option(UTOPIA_ENABLE_VC "SIMD wrapper" ON)
option(UTOPIA_ENABLE_VC "SIMD wrapper" OFF)
option(UTOPIA_STATIC_DEPENDENCIES_ONLY
       "Allows to restrict the linking to static libraries" OFF)

option(UTOPIA_ENABLE_MPI "Enable the cuda backend" ON)
# option(UTOPIA_ENABLE_OPENCL "Enable the opencl backend" OFF)
option(UTOPIA_ENABLE_PETSC "Enable the petsc backend" ON)
option(UTOPIA_ENABLE_SLEPC "Enable the slepc backend" ON)
option(UTOPIA_ENABLE_TRILINOS "Enable the Trilinos backend" ON)
option(UTOPIA_ENABLE_KOKKOS_SIMD "Enable kokkos intriniscs wrapper" OFF)
option(UTOPIA_ENABLE_BLAS "Enable the blas backend" ON)
option(UTOPIA_ENABLE_M3ELINSOL
       "Enable the m3elinsol algebraic multigrid backend" OFF)

option(UTOPIA_ENABLE_TRACE "enables utopia tracing facilities for regions" OFF)
option(UTOPIA_ENABLE_TRACE_EXPR
       "enables utopia tracing facilities for every expression" OFF)

option(UTOPIA_ENABLE_NO_ALLOC_REGIONS"enables utopia alloc-check facilities"
       OFF)
option(UTOPIA_ENABLE_SANITIZER "Enable clang -fsanitize=address flag" OFF)
# option(UTOPIA_ENABLE_PASSO_EXTENSIONS "Enable non-standard petsc solvers
# developed in the PASSO library" OFF)
option(UTOPIA_ENABLE_CXX14_FEATURES "Enable usage of cxx14 standard" ON)
option(UTOPIA_ENABLE_MOOSE_ENV_MODE
       "Allows to use the moose installation compilers" OFF)
option(UTOPIA_ENABLE_DEPRECATED_API
       "Decprecated functionality that will be removed in the future" ON)
option(UTOPIA_ENABLE_CLANG_TIDY "Use clang tidy static analyzer" OFF)
option(
    UTOPIA_ENABLE_PETSC_DM_PLEX
    "For the moment is to allow the \"make petsc\" command to install DMPlex dependencies"
    ON)
option(UTOPIA_PULL_REQUEST_MODE
       "Mode for ensuring that pull request standards are respected" OFF)

option(UTOPIA_ENABLE_EIGEN_3 "Look for eigen for comparing performance" OFF)
option(UTOPIA_BUILD_DOCUMENTATION
       "Use Doxygen to create the HTML based API documentation" ON)

option(UTOPIA_SCRIPTING
       "Enable exports for other languages (C, Python, etc...)" OFF)

option(UTOPIA_INSTALL_TRILINOS "Install trilinos directly" ON)
option(UTOPIA_INSTALL_PETSC "Install petsc directly" ON)
option(UTOPIA_INSTALL_PETSC_DEBUG "Install petsc directly" ON)
option(
    UTOPIA_ENABLE_GPERFTOOLS
    "use Google perf tools (https://gperftools.github.io/gperftools/cpuprofile.html)"
    OFF)

option(UTOPIA_SPACK_MODE "Spack env hacks" OFF)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)

# set(UTOPIA_DEFAULT_TARGETS)
include(CodeCoverage)

# ##############################################################################
# ##############################################################################
# ##############################################################################
set(UTOPIA_DEV_FLAGS
    "-Wall -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder -Werror=return-type" # -Werror=uninitialized
)

if(APPLE)
    set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=nonportable-include-path")
endif()

if(CMAKE_CXX_COMPILER MATCHES "nvcc")
    string(REPLACE "-Werror=reorder" "" UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS}")
    message(
        STATUS
            "Using nvcc so we remove the -Werror-reorder flag due to a bug in nvcc"
    )
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wno-clobbered -Wno-vla -Wno-pragmas
# -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-literal-suffix
# -Wno-deprecated-declarations -Wno-misleading-indentation
# -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-nonnull-compare
# -Wno-address -Wno-inline -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS -Werror"
# CACHE STRING "Warnings as errors settings")

# set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS}
# -Werror=inconsistent-missing-override")

if(UTOPIA_ENABLE_CLANG_TIDY)
    # FIXME
    find_program(
        CLANG_TIDY
        NAMES clang-tidy-mp-7.0 clang-tidy-7 clang-tidy-6.0 clang-tidy-5.0
              clang-tidy-4.0 clang-tidy
        PATHS /opt/local/bin)

    set(CMAKE_CXX_CLANG_TIDY
        "${CLANG_TIDY};-format-style=file;-fix;-export-fixes=fixes.yml"
    )# -warnings-as-errors=*; #-header-filter=.
endif()

# More annoying
set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Wextra ")

# More restrictive
if(UTOPIA_PULL_REQUEST_MODE)
    set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=deprecated-declarations ")
    # set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=unused-variable
    # -Werror=unused-local-typedef ")
endif()

# include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/external/GSL/include)

if(UTOPIA_ENABLE_DEPRECATED_API)
    set(UTOPIA_DEPRECATED_API ON)
endif()

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(UTOPIA_ENABLE_GPERFTOOLS)
    find_package(Gperftools)

    if(Gperftools_FOUND)
        link_libraries(gperftools::profiler)
    else()
        message(WARNING "GPERFTOOLS NOT FOUND")
    endif()
endif()

if(UTOPIA_ENABLE_TRACE_EXPR)
    set(UTOPIA_ENABLE_TRACE ON)
    set(UTOPIA_TRACE_EXPR_ENABLED ON)
    set(UTOPIA_TRACE_ENABLED ON)
endif()

if(UTOPIA_ENABLE_TRACE)
    set(UTOPIA_TRACE_ENABLED ON)
endif()

# if(UTOPIA_ENABLE_PASSO_EXTENSIONS) set(WITH_PASSO_EXTENSIONS TRUE) endif()

if(WIN32)
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG}   -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
endif()

message(
    STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if(UTOPIA_ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")

        set(UTOPIA_DEV_FLAGS
            "${UTOPIA_DEV_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls "
        )

    endif()
endif()

if(UTOPIA_ENABLE_LOCK_CHECKING)
    set(UTOPIA_ENABLE_LOCK_CHECK TRUE)
endif()

if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()
endif()

if(UTOPIA_ENABLE_CXX14_FEATURES)
    set(UTOPIA_WITH_CPP14 TRUE)
endif()

message(STATUS "[Status] UTOPIA_ROOT_PATH: ${UTOPIA_ROOT_PATH}")
include(${UTOPIA_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE)
    message(
        FATAL_ERROR
            "In-source builds are not allowed.
        Please create a directory and run cmake from there, passing the path
        to this source directory as the last argument.
        This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
        Please delete them.")
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################
# Utitlity functions

function(create_absolute_paths root paths absolutePaths)
    foreach(path ${${paths}})
        set(temp; ${temp} ${${root}}/${path})
    endforeach()
    set(${absolutePaths}
        ${temp}
        PARENT_SCOPE)
endfunction()

# ##############################################################################
macro(print_dependency_table)
    set(DEP_TABLE
        "\n____________________________________________________\n\n   BACKENDS and STATUS TABLE\n"
    )
    set(DEP_TABLE
        "${DEP_TABLE}--------------------------------------------------\n")
    set(DEP_TABLE
        "${DEP_TABLE}backend\t\t| active\t| developement status\n--------------------------------------------------\n"
    )
    set(DEP_TABLE "${DEP_TABLE}petsc\t\t| ${UTOPIA_PETSC}\t\t| usable\n")
    set(DEP_TABLE "${DEP_TABLE}slepc\t\t| ${UTOPIA_SLEPC}\t\t| usable\n")
    set(DEP_TABLE "${DEP_TABLE}blas\t\t| ${UTOPIA_BLAS}\t\t| usable\n")
    set(DEP_TABLE
        "${DEP_TABLE}trilinos\t| ${UTOPIA_TRILINOS}\t\t| basic algebra usable\n"
    )
    # SET(DEP_TABLE "${DEP_TABLE}openblas\t| ${OPEN_BLAS_FOUND}\t\t| usable\n")
    set(DEP_TABLE "${DEP_TABLE}lapack\t\t| ${UTOPIA_LAPACK}\t\t| usable\n")
    set(DEP_TABLE "${DEP_TABLE}umfpack\t\t| ${UTOPIA_UMFPACK}\t\t| disabled\n")
    # SET(DEP_TABLE "${DEP_TABLE}cuda+thrust\t| ${UTOPIA_CUDA}\t\t| not usable
    # but can be easily developed further\n") SET(DEP_TABLE
    # "${DEP_TABLE}opencl\t\t| ${UTOPIA_OPENCL}\t\t| prototyping\n")
    set(DEP_TABLE
        "${DEP_TABLE}m3elinsol\t| ${UTOPIA_M3ELINSOL}\t\t| work in progress\n")
    set(DEP_TABLE
        "${DEP_TABLE}____________________________________________________\n")
    message(STATUS ${DEP_TABLE})
endmacro()

# ##############################################################################
function(utopia_add_library libraryRootDir subDirs)
    set(THIS_MODULE_SOURCES "")
    set(THIS_MODULE_HEADERS "")
    find_project_files(${libraryRootDir} "${subDirs}" THIS_MODULE_HEADERS
                       THIS_MODULE_SOURCES)
    target_sources(utopia PRIVATE ${THIS_MODULE_SOURCES})
    install(FILES ${THIS_MODULE_HEADERS} DESTINATION include)
    foreach(MODULE ${subDirs})
        target_include_directories(
            utopia BEFORE PUBLIC $<BUILD_INTERFACE:${libraryRootDir}/${MODULE}>)
    endforeach(MODULE)
endfunction()

# ##############################################################################

function(find_project_files rootPath dirPaths headers sources)
    set(verbose TRUE)

    set(theaders ${${headers}})
    set(tsources ${${sources}})

    foreach(INCLUDE_PATH ${dirPaths})
        # INCLUDE_DIRECTORIES(${rootPath}/${INCLUDE_PATH})

        file(GLOB TEMP_HPPSRC CONFIGURE_DEPENDS
             "${rootPath}/${INCLUDE_PATH}/*.cpp")
        file(GLOB TEMP_SRC CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.c")
        file(GLOB TEMP_HPPDR CONFIGURE_DEPENDS
             "${rootPath}/${INCLUDE_PATH}/*.hpp")
        file(GLOB TEMP_HDR CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.h")

        source_group(
            ${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC};
                                  ${TEMP_SRC}; ${TEMP_UI})

        set(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
        set(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
    endforeach(INCLUDE_PATH)

    set(${headers}
        ${theaders}
        PARENT_SCOPE)
    set(${sources}
        ${tsources}
        PARENT_SCOPE)
endfunction()

# ##############################################################################

function(append_list_to_string_for_makefile OUTPUT_STRING PRE_LIST POST_LIST
         ITEM_PREFIX LIST_OF_STRINGS)
    # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro
    list(REMOVE_DUPLICATES LIST_OF_STRINGS)

    set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    foreach(LIST_EL ${LIST_OF_STRINGS})

        string(REGEX MATCH "(.*)framework" IS_FRAME_WORK "${LIST_EL}")
        if(IS_FRAME_WORK)
            get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
            string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH
                                 "${LIST_EL}")
            string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME
                                 "${FRAMEWORK_NAME}")
            set(TEMPSTR
                "${TEMPSTR}-F${FRAMEWORK_PATH} -framework ${FRAMEWORK_NAME} ")
        else()
            set(TEMPSTR "${TEMPSTR}${ITEM_PREFIX}${LIST_EL} ")
        endif()
    endforeach(LIST_EL)

    set(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    set(${OUTPUT_STRING}
        ${TEMPSTR}
        PARENT_SCOPE)
endfunction()

# ##############################################################################

function(append_list_to_string_for_cmake OUTPUT_STRING PRE_LIST POST_LIST
         LIST_OF_STRINGS)
    # http://stackoverflow.com/questions/5248749/passing-a-list-to-a-cmake-macro

    list(REMOVE_DUPLICATES LIST_OF_STRINGS)

    set(TEMPSTR ${${OUTPUT_STRING}} ${PRE_LIST})

    foreach(LIST_EL ${LIST_OF_STRINGS})
        set(TEMPSTR "${TEMPSTR}\t\"${LIST_EL}\"\n")
    endforeach(LIST_EL)

    set(TEMPSTR "${TEMPSTR} ${POST_LIST}")
    set(${OUTPUT_STRING}
        ${TEMPSTR}
        PARENT_SCOPE)
endfunction()

# ##############################################################################
# ##############################################################################

set(UTOPIA_WITH_CPP11 ON)

# ##############################################################################
# ##############################################################################
# ##############################################################################

find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
    message(
        WARNING
            "Doxygen is needed to build the documentation. Please install it correctly"
    )
endif()
# -- Configure the Template Doxyfile for our specific project
configure_file(${UTOPIA_ROOT_PATH}/Doxyfile.txt ${CMAKE_BINARY_DIR} @ONLY
               IMMEDIATE)
# -- Add a custom target to run Doxygen when ever the project is built
add_custom_target(
    docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile.txt
    SOURCES ${CMAKE_BINARY_DIR}/Doxyfile.txt)
# IF you do NOT want the documentation to be generated EVERY time you build the
# project then leave out the 'ALL' keyword from the above command.

# ##############################################################################
# ##############################################################################
# ##############################################################################

if(UTOPIA_INSTALL_TRILINOS)
    include(InstallTrilinos)
endif()

if(UTOPIA_INSTALL_PETSC)
    include(InstallPetsc)
endif()

if(UTOPIA_INSTALL_PETSC_DEBUG)
    include(InstallPetscDebug)
endif()

# using this to defining UTOPIA_LAMBDA correctly
if(UTOPIA_ENABLE_TRILINOS)
    find_package(Trilinos)
    if(Trilinos_FOUND)
        # include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS}
        # ${Trilinos_TPL_INCLUDE_DIRS}) set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
        # ${Kokkos_CXX_FLAGS}") message(STATUS
        # "Kokkos_CXX_FLAGS=${Kokkos_CXX_FLAGS},${CMAKE_CXX_FLAGS}")
        # message(STATUS "HERE:
        # ${Trilinos_CXX_COMPILER_FLAGS},${Kokkos_CXX_FLAGS},${TRILINOS_DEFINITIONS},${CMAKE_CXX_STANDARD},${INTERFACE_COMPILE_FEATURES}")
        set(UTOPIA_WITH_TRILINOS ON)
    endif()
endif()

list(APPEND UTOPIA_HEADERS utopia.hpp utopia_Base.hpp utopia_Instance.hpp
     utopia_DeprecatedHeaders.hpp)

list(APPEND UTOPIA_SOURCES utopia.cpp utopia_Instance.cpp utopia_Base.cpp)
add_library(utopia ${UTOPIA_SOURCES})

if(UTOPIA_ENABLE_KOKKOS_SIMD)
    set(UTOPIA_WITH_KOKKOS_SIMD TRUE)

    target_include_directories(
        utopia PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external/simd-math>)
endif()

# We use this to force c++14
set_utopia_compiler_features(utopia)

utopia_link_default_targets(utopia)
target_include_directories(
    utopia BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                         $<INSTALL_INTERFACE:include>)
target_include_directories(
    utopia BEFORE PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(
    TARGETS utopia
    EXPORT UtopiaTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

add_subdirectory(examples)
add_subdirectory(backend/polymorphic)

if(UTOPIA_ENABLE_VC)
    add_subdirectory(backend/vc)
endif()

# backend modules
if(UTOPIA_ENABLE_BLAS)
    add_subdirectory(backend/blas)
endif()

if(UTOPIA_ENABLE_TRILINOS)
    add_subdirectory(backend/trilinos)
endif()

if(UTOPIA_ENABLE_PETSC)
    add_subdirectory(backend/petsc)
endif()

if(UTOPIA_ENABLE_OPENCL)
    add_subdirectory(backend/opencl)
endif()

if(UTOPIA_ENABLE_M3ELINSOL)
    add_subdirectory(backend/m3elinsol)
endif()

if(UTOPIA_SCRIPTING)
    add_subdirectory(scripting)
endif()

# weird stuff goes here
if(USE_SPIKE_SOLVERS)
    list(APPEND UTOPIA_MODULES spike/solvers)
endif()

# core modules
set(UTOPIA_MODULES
    interoperability
    backend
    solvers
    differentiable
    core
    ui
    utils)

foreach(MODULE ${UTOPIA_MODULES})
    add_subdirectory(${MODULE})
endforeach(MODULE)

add_definitions(${UTOPIA_DEFS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")

# config is created based on the successful backend compilations
configure_file(${UTOPIA_ROOT_PATH}/utopia_Config.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Config.hpp)
configure_file(${UTOPIA_ROOT_PATH}/utopia_Version.hpp.in
               ${CMAKE_BINARY_DIR}/utopia_Version.hpp)

# ##############################################################################
include(WriteCompilerDetectionHeader)

write_compiler_detection_header(
    FILE utopia_compiler_features.hpp PREFIX UTOPIA
    COMPILERS GNU Clang MSVC Intel AppleClang
    FEATURES # cxx_final #already in c++11 cxx_noexcept #already in c++11
             cxx_attribute_deprecated # in c++14
    ALLOW_UNKNOWN_COMPILERS)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/utopia_compiler_features.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/utopia_Version.hpp
        DESTINATION include)

# ##############################################################################

# test suite
add_executable(utopia_test EXCLUDE_FROM_ALL test.cpp)
utopia_link_default_targets(utopia_test)
target_link_libraries(utopia_test PRIVATE utopia)
include(tests/tests.cmake)

# bench suite
add_executable(utopia_bench EXCLUDE_FROM_ALL bench.cpp)
utopia_link_default_targets(utopia_bench)
target_link_libraries(utopia_bench PRIVATE utopia)
include(benchmarks/benchmarks.cmake)

# apps suite
add_executable(utopia_exec EXCLUDE_FROM_ALL main.cpp)
utopia_link_default_targets(utopia_exec)
target_link_libraries(utopia_exec PRIVATE utopia)
include(apps/apps.cmake)

install(FILES ${UTOPIA_ROOT_PATH}/utopia.hpp DESTINATION include)
install(FILES ${UTOPIA_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/utopia_Config.hpp DESTINATION include)

# Generate files useful for dependent build systems
list(REVERSE UTOPIA_THIRDPARTY_LIBRARIES)
list(REMOVE_DUPLICATES UTOPIA_THIRDPARTY_LIBRARIES)
list(REVERSE UTOPIA_THIRDPARTY_LIBRARIES)
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_makefile(
    CONFIG_CONTENT "UTOPIA_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include " "\n"
    "-I" "${UTOPIA_THIRDPARTY_INCLUDES}")
append_list_to_string_for_makefile(
    CONFIG_CONTENT "UTOPIA_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib "
    " -lutopia" "-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_makefile(CONFIG_CONTENT " " "\n" ""
                                   "${UTOPIA_THIRDPARTY_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_DEFS=" "\n" ""
                                   "${UTOPIA_DEFS}")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_FLAGS=${CMAKE_CXX_FLAGS}\n")
set(CONFIG_CONTENT
    "${CONFIG_CONTENT}UTOPIA_CXX_COMPILER=${CMAKE_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_C_COMPILER=${CMAKE_C_COMPILER}\n")

file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(
    CONFIG_CONTENT
    "LIST(APPEND UTOPIA_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n"
    "\n)\n\n" "${UTOPIA_THIRDPARTY_INCLUDES}")
append_list_to_string_for_makefile(
    CONFIG_CONTENT
    "LIST(APPEND UTOPIA_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib\n\t-lutopia"
    "" "\n\t-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" ")\n\n"
                                "${UTOPIA_THIRDPARTY_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_DEFS\n"
                                "\n)\n\n" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_FLAGS\n"
                                "\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(
    CONFIG_CONTENT "LIST(APPEND UTOPIA_C_COMPILER\n" "\n)\n\n"
    "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(
    CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_COMPILER\n" "\n)\n\n"
    "${CMAKE_CXX_COMPILER}")
file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.cmake ${CONFIG_CONTENT})
list(JOIN CONFIG_CONTENT "\n" OLD_CMAKE_CONFIG_CONTENT)

install(FILES ${CMAKE_BINARY_DIR}/utopia-config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia-config.cmake DESTINATION config)

# new style exported targets, this allows to use target_link_libraries(my_target
# Utopia::utopia)
install(
    EXPORT UtopiaTargets
    FILE UtopiaTargets.cmake
    NAMESPACE Utopia::
    DESTINATION lib/cmake)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtopiaConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.cmake"
    INSTALL_DESTINATION "lib/cmake")
# generate the version file for the config file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfigVersion.cmake"
    VERSION "${Utopia_VERSION_MAJOR}.${Utopia_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion)
# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/UtopiaConfig.cmake
        DESTINATION lib/cmake/)

get_target_property(XX_LINK_DIRECTORIES utopia_exec LINK_DIRECTORIES)
get_target_property(XX_LINK_LIBRARIES utopia_exec LINK_INTERFACE_LIBRARIES)
get_target_property(XX_INCLUDE_DIRECTORIES utopia_exec
                    INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "LINK_DIRECTORIES=${XX_LINK_DIRECTORIES}")
message(STATUS "LINK_LIBRARIES=${XX_LINK_LIBRARIES}")
message(STATUS "INCLUDE_DIRECTORIES=${XX_INCLUDE_DIRECTORIES}")

print_dependency_table()

# find_package(OpenMP) if(OPENMP_FOUND) SET (CMAKE_Fortran_FLAGS
# "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}") SET (CMAKE_CXX_FLAGS
# "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}") SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}
# ${OpenMP_C_FLAGS}") set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS}
# ${OpenMP_CXX_FLAGS}") endif()

# ##############################################################################

# enable_testing()
add_custom_target(utest COMMAND time ${CMAKE_BINARY_DIR}/utopia_test -verbose)
# add_test(all_tests ${CMAKE_BINARY_DIR}/utopia_test -verbose)
add_dependencies(utest utopia_test)

# ##############################################################################

add_custom_target(bench COMMAND time ${CMAKE_BINARY_DIR}/utopia_bench -verbose)
add_dependencies(bench utopia_bench)

# ##############################################################################

add_custom_target(bench_and_test)
add_dependencies(bench_and_test bench utest)

# ##############################################################################

add_custom_target(complete)
add_dependencies(complete utopia utopia_test utopia_bench utopia_exec
                 utopia_examples)

# ##############################################################################
target_format_setup(utopia_test)
# ##############################################################################
