cmake_minimum_required(VERSION 2.8)
project(utopia_fe)

#To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)


#################################################################################
#################################################################################
#################################################################################

set(UTOPIA_FE_ROOT_PATH ${CMAKE_SOURCE_DIR})
set(UTOPIA_FE_SRC_PATH ${UTOPIA_FE_ROOT_PATH}/src)
set(UTOPIA_FE_CMAKES_PATH ${UTOPIA_FE_ROOT_PATH}/libcmake)
set(TYPE DYNAMIC)

if (ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
		SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ")
	endif()
endif()

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
################################## Include cmake functions and find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_FE_ROOT_PATH}/cmake")
file(GLOB_RECURSE FUNC_LIB "{UTOPIA_FE_ROOT_PATH}/cmake/utils" )
foreach(funcLib ${FUNC_LIB})
	include(${funcLib})
endforeach(funcLib)

find_package(Utopia)
if(UTOPIA_FOUND)
    list(APPEND CMAKE_CXX_FLAGS "-std=c++11")
    include_directories(${UTOPIA_INCLUDES})
    add_definitions(${UTOPIA_DEFS})
else()
    message(FATAL_ERROR "Utopia not found.")
endif()

find_package(LIBMESH)
if(LIBMESH_FOUND)
	include_directories(${LIBMESH_INCLUDE_DIR})
	list(APPEND UTOPIA_FE_LIBRARIES ${LIBMESH_LIBRARIES})
	# list(APPEND FELANG_LINK_DIRS${POSTFIX} ${PARMETIS_LIB_DIR})
else()
	message(FATAL_ERROR "could not find LibMesh")
endif()

find_package(Boost)
find_package(BoostExtended)
if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS})
		link_directories(${Boost_LIBRARY_DIRS})
		list(APPEND UTOPIA_FE_LIBRARIES ${Boost_LIBRARIES})
else()
	MESSAGE(WARNING "Boost not found.")
endif()

add_definitions("-D_LIBCPP_DEBUG=0")

list(APPEND HEADERS utopia_fe.hpp)
list(APPEND HEADERS rapidxml/rapidxml.hpp rapidxml/rapidxml_iterators.hpp rapidxml/rapidxml_print.hpp rapidxml/rapidxml_utils.hpp)
set(SOURCES "")


include_directories(rapidxml)
include_directories(transfer)

list(APPEND SOURCES utopia_FEDSLBaseExamples.cpp)
list(APPEND SOURCES utopia_FEDSLMortarExamples.cpp)
list(APPEND SOURCES utopia_FEDSLLeastSquaresExamples.cpp)
list(APPEND SOURCES utopia_MixedFESpaceExample.cpp)
list(APPEND SOURCES utopia_TimeDiffExamples.cpp)
list(APPEND SOURCES utopia_Biomechanics.cpp)
list(APPEND SOURCES utopia_UGMeshReader.cpp)
list(APPEND SOURCES utopia_NormalTangentialCoordinateSystem.cpp)
list(APPEND SOURCES utopia_GeometryTest.cpp)

add_subdirectory(transfer)
add_library(utopia_fe STATIC ${SOURCES} ${HEADERS})
target_link_libraries(utopia_fe ${UTOPIA_LIBRARIES} ${UTOPIA_FE_LIBRARIES} ${TRANSFER_LIB})

add_executable(utopia_fe_exec  main.cpp ${HEADERS})
target_link_libraries(utopia_fe_exec utopia_fe)

