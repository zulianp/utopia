cmake_minimum_required(VERSION 3.5)

project(utopia)
# cmake_policy(SET "CMP0079" NEW)

set(UTOPIA_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/cmake")

#################################################################################
######################### UTOPIA VERSION INFORMATION ############################
#################################################################################

include(UtopiaGitVersion)
set(UTOPIA_VERSION_MAJOR 0)
set(UTOPIA_VERSION_MINOR 1)
set(UTOPIA_VERSION_PATCH 0)
set(UTOPIA_VERSION "${UTOPIA_VERSION_MAJOR}.${UTOPIA_VERSION_MINOR}.${UTOPIA_VERSION_PATCH}")
#################################################################################





if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

#To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)

#################################################################################
#################################################################################
#################################################################################

option(UTOPIA_ARCHIVE_ONLY "Allows to compile only the archive file. Suitable for linux" OFF)
option(UTOPIA_STATIC_DEPENDENCIES_ONLY "Allows to restrict the linking to static libraries" OFF)

option(TRY_WITH_MPI "Enable the cuda backend" ON)
# option(TRY_WITH_OPENCL "Enable the opencl backend" OFF)
option(TRY_WITH_PETSC "Enable the petsc backend" ON)
option(TRY_WITH_SLEPC "Enable the slepc backend" ON)
option(TRY_WITH_TRILINOS "Enable the Trilinos backend" ON)
option(TRY_WITH_BLAS "Enable the blas backend" ON)
option(TRY_WITH_M3ELINSOL "Enable the m3elinsol algebraic multigrid backend" OFF)

option(ENABLE_TRACE "enables utopia tracing facilities" OFF)
option(ENABLE_SANITIZER "Enable clang -fsanitize=address flag" OFF)
option(ENABLE_PASSO_EXTENSIONS "Enable non-standard petsc solvers developed in the PASSO library" OFF)
option(ENABLE_CXX14_FEATURES "Enable usage of cxx14 standard" OFF)
option(ENABLE_MOOSE_ENV_MODE "Allows to use the moose installation compilers" ON)

option(TRY_WITH_EIGEN_3 "Look for eigen for comparing performance" OFF)
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)

#################################################################################
#################################################################################
#################################################################################
set(UTOPIA_DEV_FLAGS
    "-Wall -Werror=enum-compare -Werror=uninitialized -Werror=delete-non-virtual-dtor -Werror=reorder "

)

# set(CMAKE_CXX_FLAGS "-Wall -Wno-clobbered -Wno-vla -Wno-pragmas -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-deprecated-declarations -Wno-misleading-indentation -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-nonnull-compare -Wno-address -Wno-inline -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS -Werror" CACHE STRING "Warnings as errors settings")

# set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=inconsistent-missing-override")

# More annoying
set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Wextra ")

# More restrictive
# set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=unused-variable -Werror=unused-local-typedef -Werror=deprecated-declarations ")

include_directories(${CMAKE_SOURCE_DIR}/external/GSL/include)


if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()



if(ENABLE_TRACE)
    set(UTOPIA_TRACE_ENABLED ON)
endif()

if(ENABLE_PASSO_EXTENSIONS)
    set(WITH_PASSO_EXTENSIONS TRUE)
endif()

IF(WIN32)
    SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
ENDIF()

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if (ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")

        SET(UTOPIA_DEV_FLAGS  "${UTOPIA_DEV_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ")

    endif()
endif()

if(ENABLE_LOCK_CHECKING)
    set(ENABLE_LOCK_CHECK TRUE)
endif()

if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
    IF(WIN32)
         SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    ELSE()
         SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    ENDIF()
endif()

if(ENABLE_CXX14_FEATURES)
    SET(WITH_CPP14 TRUE)
endif()



message(STATUS "[Status] UTOPIA_ROOT_PATH: ${UTOPIA_ROOT_PATH}")
include(${UTOPIA_ROOT_PATH}/cmake/UtopiaCompilerFeatures.cmake)

if( CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE )
    message(FATAL_ERROR
        "In-source builds are not allowed.
        CMake would overwrite the makefiles distributed with Moose.
        Please create a directory and run cmake from there, passing the path
        to this source directory as the last argument.
        This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
        Please delete them."
    )
endif()

#################################################################################
#################################################################################
#################################################################################
# Utitlity functions

function(create_absolute_paths root paths absolutePaths)
    FOREACH(path ${${paths}})
    SET(temp; ${temp} ${${root}}/${path})
    ENDFOREACH()
    SET(${absolutePaths} ${temp} PARENT_SCOPE)
endfunction()

#################################################################################
macro(print_dependency_table)
    SET(DEP_TABLE "\n____________________________________________________\n\n   BACKENDS and STATUS TABLE\n")
    SET(DEP_TABLE "${DEP_TABLE}--------------------------------------------------\n")
    SET(DEP_TABLE "${DEP_TABLE}backend\t\t| active\t| developement status\n--------------------------------------------------\n")
    SET(DEP_TABLE "${DEP_TABLE}petsc\t\t| ${UTOPIA_PETSC}\t\t| usable\n")
    SET(DEP_TABLE "${DEP_TABLE}slepc\t\t| ${UTOPIA_SLEPC}\t\t| usable\n")
    SET(DEP_TABLE "${DEP_TABLE}blas\t\t| ${UTOPIA_BLAS}\t\t| usable\n")
    SET(DEP_TABLE "${DEP_TABLE}trilinos\t| ${UTOPIA_TRILINOS}\t\t| basic algebra usable\n")
    # SET(DEP_TABLE "${DEP_TABLE}openblas\t| ${OPEN_BLAS_FOUND}\t\t| usable\n")
    SET(DEP_TABLE "${DEP_TABLE}lapack\t\t| ${UTOPIA_LAPACK}\t\t| usable\n")
    SET(DEP_TABLE "${DEP_TABLE}umfpack\t\t| ${UTOPIA_UMFPACK}\t\t| disabled\n")
    # SET(DEP_TABLE "${DEP_TABLE}cuda+thrust\t| ${UTOPIA_CUDA}\t\t| not usable but can be easily developed further\n")
    # SET(DEP_TABLE "${DEP_TABLE}opencl\t\t| ${UTOPIA_OPENCL}\t\t| prototyping\n")
    SET(DEP_TABLE "${DEP_TABLE}m3elinsol\t| ${UTOPIA_M3ELINSOL}\t\t| work in progress\n")
    SET(DEP_TABLE "${DEP_TABLE}____________________________________________________\n")
    message(STATUS ${DEP_TABLE})
endmacro()

#################################################################################

function(find_project_files rootPath dirPaths headers sources)
    SET(verbose TRUE)

    SET(theaders ${${headers}})
    SET(tsources ${${sources}})

    SET(ROOT_PATH ${${rootPath}})
    SET(DIR_PATHS ${${dirPaths}})

    FOREACH(INCLUDE_PATH ${DIR_PATHS})
        INCLUDE_DIRECTORIES(${ROOT_PATH}/${INCLUDE_PATH})

        file(GLOB TEMP_HPPSRC "${ROOT_PATH}/${INCLUDE_PATH}/*.cpp" )
        file(GLOB TEMP_SRC "${ROOT_PATH}/${INCLUDE_PATH}/*.c" )
        file(GLOB TEMP_HPPDR "${ROOT_PATH}/${INCLUDE_PATH}/*.hpp" )
        file(GLOB TEMP_HDR "${ROOT_PATH}/${INCLUDE_PATH}/*.h" )

        SOURCE_GROUP(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})

        SET(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
        SET(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
    ENDFOREACH(INCLUDE_PATH)

    SET(${headers} ${theaders} PARENT_SCOPE)
    SET(${sources} ${tsources} PARENT_SCOPE)
endfunction()

#################################################################################

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

#################################################################################
#################################################################################

set(WITH_CPP11 ON)


#################################################################################
#################################################################################
#################################################################################

find_package(Doxygen)
IF (NOT DOXYGEN_FOUND)
message(WARNING "Doxygen is needed to build the documentation. Please install it correctly")
ENDIF()
#-- Configure the Template Doxyfile for our specific project
configure_file(${UTOPIA_ROOT_PATH}/Doxyfile.txt  ${CMAKE_BINARY_DIR}  @ONLY IMMEDIATE)
#-- Add a custom target to run Doxygen when ever the project is built
add_custom_target(docs COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile.txt SOURCES ${CMAKE_BINARY_DIR}/Doxyfile.txt)
# IF you do NOT want the documentation to be generated EVERY time you build the project then leave out the 'ALL' keyword from the above command.

#################################################################################
#################################################################################
#################################################################################

# using this to defining UTOPIA_LAMBDA correctly
if(TRY_WITH_TRILINOS)
    find_package(Trilinos)
    if(Trilinos_FOUND)
        include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
        set(WITH_TRILINOS ON)
    endif()
endif()

add_subdirectory(backend/polymorphic)

# backend modules
if(TRY_WITH_BLAS)
    add_subdirectory(backend/blas)
endif()

if(TRY_WITH_TRILINOS)
    add_subdirectory(backend/trilinos)
endif()

if(TRY_WITH_PETSC)
    add_subdirectory(backend/petsc)
endif()

if(TRY_WITH_OPENCL)
    add_subdirectory(backend/opencl)
endif()

if(TRY_WITH_M3ELINSOL)
    add_subdirectory(backend/m3elinsol)
endif()

# weird stuff goes here
if(USE_SPIKE_SOLVERS)
    list(APPEND UTOPIA_MODULES spike/solvers)
endif()


# core modules
set(UTOPIA_MODULES interoperability backend solvers differentiable core ui utils)

foreach(MODULE ${UTOPIA_MODULES})
    add_subdirectory(${MODULE})
endforeach(MODULE)



add_definitions(${UTOPIA_DEFS})
include_directories(${UTOPIA_INCLUDES})
include_directories(${UTOPIA_INCLUDE_PATHS})
include_directories(${UTOPIA_MODULES})
include_directories(${CMAKE_BINARY_DIR})
include_directories(.)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_DEV_FLAGS}")

list(APPEND UTOPIA_HEADERS
    utopia.hpp
    utopia_Base.hpp
    utopia_Instance.hpp
    )

list(APPEND UTOPIA_SOURCES
    utopia.cpp
    utopia_Instance.cpp
    utopia_Base.cpp
    )

# config is created based on the successful backend compilations
configure_file(${UTOPIA_ROOT_PATH}/utopia_Config.hpp.in ${CMAKE_BINARY_DIR}/utopia_Config.hpp)
configure_file(${UTOPIA_ROOT_PATH}/utopia_Version.hpp.in ${CMAKE_BINARY_DIR}/utopia_Version.hpp)


#########################################################
include(WriteCompilerDetectionHeader)

write_compiler_detection_header(
  FILE utopia_compiler_features.hpp
  PREFIX UTOPIA
  COMPILERS GNU Clang MSVC Intel AppleClang
  FEATURES
             # cxx_final #already in c++11
             # cxx_noexcept #already in c++11
             cxx_attribute_deprecated # in c++14
  ALLOW_UNKNOWN_COMPILERS
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/utopia_compiler_features.hpp
    DESTINATION include
)

##########################################################


add_library(utopia STATIC utopia.cpp ${UTOPIA_SOURCES})
target_link_libraries(utopia ${UTOPIA_COMPONENTS})
target_link_libraries(utopia ${UTOPIA_LIBRARIES})

if(NOT UTOPIA_ARCHIVE_ONLY)
    # add_library(utopia_dynamic SHARED utopia.cpp ${UTOPIA_SOURCES})
    # target_link_libraries(utopia_dynamic ${UTOPIA_COMPONENTS})
    # target_link_libraries(utopia_dynamic ${UTOPIA_LIBRARIES})

    # test suite
    add_executable(utopia_test test.cpp)
    target_link_libraries(utopia_test utopia)
    # add_subdirectory(tests)
    include(tests/tests.cmake)

    add_executable(utopia_exec main.cpp ${UTOPIA_HEADERS} ${UTOPIA_SOURCES})
    target_link_libraries(utopia_exec utopia)
    target_link_libraries(utopia_exec ${UTOPIA_COMPONENTS})
    target_link_libraries(utopia_exec ${UTOPIA_LIBRARIES})

    # set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
    
    install(TARGETS utopia utopia_exec
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
else()
    install(TARGETS utopia
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()

install(FILES ${UTOPIA_ROOT_PATH}/utopia.hpp DESTINATION include)
install(FILES ${UTOPIA_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/utopia_Config.hpp DESTINATION include)

# Generate files useful for dependent build systems
LIST(REVERSE UTOPIA_LIBRARIES)
LIST(REMOVE_DUPLICATES UTOPIA_LIBRARIES)
LIST(REVERSE UTOPIA_LIBRARIES)
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_INCLUDES=-I${CMAKE_INSTALL_PREFIX}/include " "\n" "-I" "${UTOPIA_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_LIBRARIES=-L${CMAKE_INSTALL_PREFIX}/lib " " -lutopia" "-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_makefile(CONFIG_CONTENT " " "\n" "" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_DEFS=" "\n" "" "${UTOPIA_DEFS}")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_FLAGS=${CMAKE_CXX_FLAGS}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_COMPILER=${CMAKE_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_C_COMPILER=${CMAKE_C_COMPILER}\n")


file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_INCLUDES\n\t${CMAKE_INSTALL_PREFIX}/include\n" 		"\n)\n\n" "${UTOPIA_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "LIST(APPEND UTOPIA_LIBRARIES\n\t-L${CMAKE_INSTALL_PREFIX}/lib" "\n\t-lutopia" "\n\t-l" "${UTOPIA_COMPONENTS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "\n" ")\n\n" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_DEFS\n" 			"\n)\n\n" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_FLAGS\n"  	"\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_C_COMPILER\n"  	"\n)\n\n" "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_COMPILER\n" 	"\n)\n\n" "${CMAKE_CXX_COMPILER}")

file(WRITE ${CMAKE_BINARY_DIR}/utopia-config.cmake ${CONFIG_CONTENT})

install(FILES ${CMAKE_BINARY_DIR}/utopia-config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia-config.cmake 	 DESTINATION config)

print_dependency_table()


# find_package(OpenMP)
# if(OPENMP_FOUND)
# 	SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
# 	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# 	SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# 	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()




#############################
### petsc-cuda stuff

link_directories(${CMAKE_CURRENT_BINARY_DIR})
