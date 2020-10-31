set(LIBMESH_INCLUDES_MY_LOCATIONS ~/Desktop/libmesh/include)
set(LIBMESH_LIB_MY_LOCATIONS ~/Desktop/libmesh/lib)
set(LIBMESH_BIN_MY_LOCATIONS ~/Desktop/libmesh/bin)

# list(APPEND LIBMESH_INCLUDES_MY_LOCATIONS
# /Users/mariagiuseppina/project_with_lib_and_moonolith/libmesh_build/include)
# list(APPEND LIBMESH_LIB_MY_LOCATIONS
# /Users/mariagiuseppina/project_with_lib_and_moonolith/libmesh_build/lib)
# list(APPEND LIBMESH_BIN_MY_LOCATIONS
# /Users/mariagiuseppina/project_with_lib_and_moonolith/libmesh_build)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE NONE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
if(${BUILD_TYPE_UPPER} MATCHES DEBUG)
    set(METHOD dbg)
else()
    set(METHOD opt)
endif()
message(STATUS "linking against ${METHOD} libmesh library")

find_package(PkgConfig QUIET)

# set(LIBMESH_DIR LIBMESH_DIR-NOTFOUND CACHE PATH "Libmesh installation
# directory")
if($ENV{LIBMESH_DIR})
    set(ENV{PKG_CONFIG_PATH}
        "$ENV{LIBMESH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

if(LIBMESH_DIR)
    set(ENV{PKG_CONFIG_PATH}
        "${LIBMESH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

pkg_check_modules(PC_LIBMESH QUIET libmesh-${METHOD})
# message(STATUS "PC_LIBMESH_FOUND = ${PC_LIBMESH_FOUND}") message(STATUS
# "PC_LIBMESH_LIBRARIES = ${PC_LIBMESH_LIBRARIES}") message(STATUS
# "PC_LIBMESH_LIBRARY_DIRS = ${PC_LIBMESH_LIBRARY_DIRS}")

# strip flags that are not definitions (-D...) message(STATUS
# "PC_LIBMESH_CFLAGS_OTHER = ${PC_LIBMESH_CFLAGS_OTHER}")
foreach(FLAG ${PC_LIBMESH_CFLAGS_OTHER})
    # message(${FLAG})
    if(${FLAG} MATCHES "^[-][D].+")
        list(APPEND PC_LIBMESH_CFLAGS_STRIPPED ${FLAG})
    endif()
endforeach()
# message(STATUS "PC_LIBMESH_CFLAGS_STRIPPED = ${PC_LIBMESH_CFLAGS_STRIPPED}")
set(LIBMESH_DEFINITIONS ${PC_LIBMESH_CFLAGS_STRIPPED})

find_path(
    LIBMESH_INCLUDE_DIR libmesh/libmesh.h
    HINTS ${LIBMESH_DIR}/include $ENV{LIBMESH_DIR}/include
          ${PC_LIBMESH_INCLUDEDIR} ${PC_LIBMESH_INCLUDE_DIRS}
          ${LIBMESH_INCLUDES_MY_LOCATIONS}
    PATH_SUFFIXES libmesh installed)

find_library(
    LIBMESH_LIBRARY
    NAMES mesh_${METHOD} mesh mesh_opt
    HINTS ${LIBMESH_DIR}/lib $ENV{LIBMESH_DIR}/lib $ENV{LIBMESH_DIR}
          ${PC_LIBMESH_LIBDIR} ${PC_LIBMESH_LIBARY_DIRS}
          ${LIBMESH_LIB_MY_LOCATIONS})

message(
    STATUS
        "HERE: $ENV{LIBMESH_DIR}/lib ${LIBMESH_LIBRARY} -> method: ${METHOD}!")

set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARY})
set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIR})

get_filename_component(LIBMESH_LIB_DIR ${LIBMESH_LIBRARY} DIRECTORY)
find_library(
    TIMPI_LIBRARY
    NAMES timpi_${METHOD} timpi timpi_opt
    HINTS ${LIBMESH_LIB_DIR})

if(TIMPI_LIBRARY)
    message(STATUS "TIMPI_LIBRARY: ${TIMPI_LIBRARY}")
    list(APPEND LIBMESH_LIBRARIES "${TIMPI_LIBRARY}")
endif()

find_program(
    LIBMESH_CONFIG_EXECUTABLE
    NAMES libmesh-config
    HINTS ${LIBMESH_DIR} ${LIBMESH_INCLUDE_DIR}/../ ${LIBMESH_BIN_MY_LOCATIONS}
          $ENV{LIBMESH_DIR}
    PATH_SUFFIXES bin
    DOC "libmesh-config executable")
mark_as_advanced(LIBMESH_CONFIG_EXECUTABLE)

exec_program(
    ${LIBMESH_CONFIG_EXECUTABLE} ARGS
    --include OUTPUT_VARIABLE
    LMC_INC_FLAG RETURN_VALUE
    LMC_INC_RET)
string(REPLACE " " ";" LMC_INC_LIST ${LMC_INC_FLAG})
foreach(IPATH ${LMC_INC_LIST})

    string(REGEX MATCH "^-I" IS_INCLUDE ${IPATH})

    if(IS_INCLUDE)
        # libmesh (~ 1.6.0-pre) added to the include list macro definitions and
        # -pthread here we filter them out
        string(REGEX REPLACE "^-I" "" IPATH ${IPATH})
        string(REGEX REPLACE "//" "/" IPATH ${IPATH})
        list(APPEND LM_INC ${IPATH})
    endif()

endforeach()
set(LIBMESH_INCLUDE_DIRS ${LM_INC})

if(PC_LIBMESH_VERSION)
    set(LIBMESH_VERSION_STRING ${PC_LIBMESH_VERSION})
endif()

# exec_program( ${LIBMESH_CONFIG_EXECUTABLE} ARGS --libs OUTPUT_VARIABLE
# LMC_LIB_FLAG RETURN_VALUE LMC_LIB_RET )

# string(REPLACE " " ";" LMC_LIB_LIST ${LMC_LIB_FLAG}) set(LIBMESH_LIBRARIES
# "${LIBMESH_LIBRARIES} ${LMC_LIB_FLAG}")

message(STATUS "LIBMESH_LIBRARIES: ${LIBMESH_LIBRARIES}")

# handle the QUIETLY and REQUIRED arguments and set LIBMESH_FOUND to TRUE if all
# listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    LIBMESH
    REQUIRED_VARS LIBMESH_LIBRARIES LIBMESH_INCLUDE_DIR
    VERSION_VAR LIBMESH_VERSION_STRING)

mark_as_advanced(LIBMESH_INCLUDE_DIR LIBMESH_LIBRARIES
                 LIBMESH_CONFIG_EXECUTABLE)
