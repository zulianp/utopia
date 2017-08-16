cmake_minimum_required(VERSION 2.8.8)
project(utopia)

#To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)
# set(CMAKE_INSTALL_RPATH "$ORIGIN")
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,$ORIGIN")
# link_directories(${CMAKE_CURRENT_BINARY_DIR})

#################################################################################
#################################################################################
#################################################################################
option(UTOPIA_ARCHIVE_ONLY "Allows to compile only the archive file. Suitable for linux" OFF)
option(UTOPIA_STATIC_DEPENDENCIES_ONLY "Allows to restrict the linking to static libraries" OFF)
option(TRY_WITH_CUDA "Allows to enable to try enabling the cuda backend" OFF)
option(TRY_WITH_OPENCL "Allows to enable to try enabling the opencl backend" OFF)
option(TRY_WITH_PETSC "Allows to enable to try enabling the petsc backend" ON)
option(TRY_WITH_BLAS "Allows to enable to try enabling the blas backend" ON)
option(TRY_WITH_EIGEN_3 "Look for eigen for comparing performance" ON)
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)
option(UTOPIA_LOG "enables utopia logging facilities" OFF)

#################################################################################
#################################################################################
#################################################################################

IF(WIN32)
	SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
	SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
ENDIF()

# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
# 	if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
# 		SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ")
# 	endif()
# endif()

if(UTOPIA_STATIC_DEPENDENCIES_ONLY)
	IF(WIN32)
	 	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
	ELSE()
	 	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	ENDIF()
endif()

set(UTOPIA_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/../cmakemodules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/cmake")
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_ROOT_PATH}/../cmakeutils")

message(STATUS "[Status] UTOPIA_ROOT_PATH: ${UTOPIA_ROOT_PATH}")

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
	SET(DEP_TABLE "${DEP_TABLE}petsc\t\t| ${PETSC_FOUND}\t\t| usable\n")
	SET(DEP_TABLE "${DEP_TABLE}blas\t\t| ${BLAS_FOUND}\t\t| usable\n")
	SET(DEP_TABLE "${DEP_TABLE}openblas\t| ${OPEN_BLAS_FOUND}\t\t| usable\n")
	SET(DEP_TABLE "${DEP_TABLE}lapack\t\t| ${LAPACK_FOUND}\t\t| usable\n")
	SET(DEP_TABLE "${DEP_TABLE}umfpack\t\t| ${UMFPACK_FOUND}\t\t| usable\n")
	SET(DEP_TABLE "${DEP_TABLE}cuda+thrust\t| ${THRUST_FOUND}\t\t| not usable but can be easily developed further\n")
	SET(DEP_TABLE "${DEP_TABLE}opencl\t\t| ${OPENCL_FOUND}\t\t| prototyping\n")
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
			# message(STATUS "here: ${LIST_EL} = ${IS_FRAME_WORK}")
			get_filename_component(FRAMEWORK_NAME "${LIST_EL}" NAME)
			string(REGEX REPLACE "/${FRAMEWORK_NAME}" "" FRAMEWORK_PATH "${LIST_EL}")
			string(REGEX REPLACE ".framework" "" FRAMEWORK_NAME "${FRAMEWORK_NAME}")

			# message(STATUS "name: ${FRAMEWORK_NAME}")
			# message(STATUS "path: ${FRAMEWORK_PATH}")

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
#################################################################################
# find all dependencies
if(TRY_WITH_BLAS)
	if(PITZ_DORA)
		find_package(OpenBLAS)
		if(OPEN_BLAS_FOUND)
			list(APPEND UTOPIA_INCLUDES ${BLAS_INCLUDE_DIR})
			list(APPEND UTOPIA_LIBRARIES ${BLAS_LIBRARIES})
			list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
			set(WITH_OPEN_BLAS ON)
		else()
			message(WARNING "[Warning] blas not found")
			set(WITH_BLAS OFF)
		endif()
	else()
		SET(OPEN_BLAS_FOUND FALSE)
		find_package(BLAS)
		if(BLAS_FOUND)
			list(APPEND UTOPIA_INCLUDES ${BLAS_INCLUDE_DIR})
			list(APPEND UTOPIA_LIBRARIES ${BLAS_LIBRARIES})
			list(APPEND UTOPIA_DEFS ${BLAS_DEFINITIONS})
			set(WITH_BLAS ON)
		else()
			SET(BLAS_FOUND FALSE)
		endif()
	endif()
endif()

#################################################################################

if(TRY_WITH_PETSC)
	
	find_package(MPI)
	if(MPI_FOUND)
		# MESSAGE( STATUS "MPI_CXX_LIBRARIES: " ${MPI_CXX_LIBRARIES} )
		# MESSAGE( STATUS "MPI_C_LIBRARIES: " ${MPI_C_LIBRARIES} )
		# MESSAGE( STATUS "MPI_FORT_LIBRARIES: " ${MPI_FORT_LIBRARIES} )
		list(APPEND UTOPIA_INCLUDES ${MPI_CXX_INCLUDE_PATH})
		list(APPEND UTOPIA_LIBRARIES  ${MPI_CXX_LIBRARIES} )
	endif()

	find_package(PETSc)
	
	if(PETSC_FOUND)
		list(APPEND UTOPIA_INCLUDES ${PETSC_INCLUDES})
		list(APPEND UTOPIA_LIBRARIES ${PETSC_LIBRARIES})

	    set(CMAKE_C_COMPILER ${PETSC_COMPILER})

	    if(NOT MPI_CXX_COMPILER)
			set(MPI_CXX_COMPILER $ENV{MPI_CXX_COMPILER})
			message(STATUS "compiler ${MPI_CXX_COMPILER}")
	    endif()

		if(MPI_CXX_COMPILER)
			set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
			set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
		else()
			 execute_process(COMMAND mpicxx -v RESULT_VARIABLE MPICXX_FAILED)

			if(MPICXX_FAILED)
				message(STATUS "Using CMAKE compiler, you can define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
			else()
	    		message(STATUS "-----------------------------------------------")
		    	message(STATUS "\n[MPI] using mpicxx for compiling c++ files.\nIf you want to use your own compiler define MPI_CXX_COMPILER=<alias_or_path_to_your_compiler>")
		    	message(STATUS "-----------------------------------------------")
		    	set(CMAKE_CXX_COMPILER mpicxx)
		    	set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
			endif()
		endif()

	    set(WITH_PETSC ON)
	    set(WITH_MPI ON)
	else()
	    message(WARNING "[Warning] PETSc not found")
	    set(WITH_PETSC OFF)
	endif()
endif()



#################################################################################
if(TRY_WITH_BLAS)
	find_package(LAPACK)
	if(LAPACK_FOUND)
	    list(APPEND UTOPIA_INCLUDES ${LAPACK_INCLUDE_DIR})
	    list(APPEND UTOPIA_LIBRARIES ${LAPACK_LIBRARIES})
	    list(APPEND UTOPIA_DEFS ${LAPACK_DEFINITIONS})
	    set(WITH_LAPACK ON)
	else()
	    message(WARNING "[Warning] lapack not found")
	    set(WITH_LAPACK OFF)
	endif()

	find_package(Umfpack)

	IF(UMFPACK_FOUND)
		list(APPEND UTOPIA_INCLUDES ${UMFPACK_INCLUDES})
		list(APPEND UTOPIA_LIBRARIES ${UMFPACK_LIBRARIES})
		set(WITH_UMFPACK ON)
	ELSE()
	    message(WARNING "[Warning] Umfpack not found")
	    set(WITH_UMFPACK OFF)
	ENDIF()
endif()

set(WITH_CPP11 ON)

#################################################################################


if(TRY_WITH_CUDA)
	find_package(CUDA)
	if(CUDA_FOUND)
		find_package(Thrust)
		if(THRUST_FOUND)

			file(GLOB CUDA_FILES "${UTOPIA_ROOT_PATH}/backend/cuda/*.cu") 
			list(APPEND CUDA_FILES backend/cuda/utopia_CUDAKernels.cu backend/cuda/utopia_Thrust.cu)
			set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER}) 
			cuda_include_directories( ${UTOPIA_ROOT_PATH}/backend/cuda ${THRUST_INCLUDE_DIR} ${UTOPIA_ROOT_PATH}/core ${UTOPIA_ROOT_PATH}/core/expressions ${CMAKE_BINARY_DIR})
			set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; "-std=c++11")
			cuda_compile(CUDA_OBJECTS ${CUDA_FILES})
			message(STATUS "${CUDA_OBJECTS}")
			cuda_add_library(utopia_cuda ${CUDA_FILES})
			set(WITH_CUDA ON)

			# list(APPEND UTOPIA_INCLUDES ${THRUST_INCLUDE_DIR})
			# message(FATAL_ERROR "${THRUST_INCLUDE_DIR}")
		else()
			message(FATAL_ERROR "Thrust not found")	
		endif()
	else()
		message(FATAL_ERROR "Cuda not found")
	endif()
else()
	set(WITH_CUDA FALSE)
endif()

if(TRY_WITH_OPENCL)
	find_package(OpenCL)
	if(OPENCL_FOUND)
		set(OPENCL_FOUND TRUE)
		set(WITH_UTOPIA_OPENCL TRUE)
		list(APPEND UTOPIA_LIBRARIES ${OPENCL_LIBRARIES})
		list(APPEND UTOPIA_INCLUDES ${OPENCL_INCLUDE_DIRS})
	else()
		set(OPENCL_FOUND FALSE)
	endif()
endif()


if(TRY_WITH_EIGEN_3)
	find_package(Eigen3)
	if(EIGEN3_FOUND)
		set(WITH_EIGEN_3 TRUE)
		list(APPEND UTOPIA_INCLUDES ${EIGEN3_INCLUDE_DIR})
	endif()
endif()	

#################################################################################
#################################################################################
#################################################################################

find_package(Doxygen)
IF (NOT DOXYGEN_FOUND)
message(WARNING "Doxygen is needed to build the documentation. Please install it correctly")
ENDIF()
#-- Configure the Template Doxyfile for our specific project
configure_file(${UTOPIA_ROOT_PATH}/Doxyfile.txt  ${UTOPIA_ROOT_PATH}/bin  @ONLY IMMEDIATE)
#-- Add a custom target to run Doxygen when ever the project is built
add_custom_target(docs COMMAND ${DOXYGEN_EXECUTABLE} ${UTOPIA_ROOT_PATH}/bin/Doxyfile.txt SOURCES ${UTOPIA_ROOT_PATH}/bin/Doxyfile.txt)
# IF you do NOT want the documentation to be generated EVERY time you build the project then leave out the 'ALL' keyword from the above command.

#################################################################################
#################################################################################
#################################################################################

set(UTOPIA_MODULES core core/expressions core/expressions_parallel core/eval backend solvers solvers/linear solvers/linear/homemade solvers/nonlinear solvers/nonlinear/hessian_approximations solvers/multilevel solvers/domain_decomposition solvers/nonlinear/line_search solvers/nonlinear/trust_region solvers/nonlinear/constrained solvers/nonlinear/constrained/quadratic_programming solvers/eig utils differentiable tests)

if(BLAS_FOUND OR OPEN_BLAS_FOUND)
	list(APPEND UTOPIA_MODULES backend/blas)
endif()

if(LAPACK_FOUND)
	list(APPEND UTOPIA_MODULES solvers/linear/lapack solvers/eig/lapack) 
endif()

if(PETSC_FOUND)
	list(APPEND UTOPIA_MODULES backend/petsc solvers/linear/petsc core/eval/petsc)
endif()

if(UMFPACK_FOUND)
	list(APPEND UTOPIA_MODULES solvers/linear/umfpack)
endif()

if(THRUST_FOUND)
	list(APPEND UTOPIA_MODULES backend/cuda)
endif()

if(OPENCL_FOUND)
	list(APPEND UTOPIA_MODULES backend/opencl)
endif()

if(USE_SPIKE_SOLVERS)
	list(APPEND UTOPIA_MODULES spike/solvers)
endif()

if(UTOPIA_LOG)
	add_definitions("-DUTOPIA_LOG_ENABLED")
endif()


create_absolute_paths(UTOPIA_ROOT_PATH UTOPIA_MODULES UTOPIA_INCLUDE_PATHS)
find_project_files(UTOPIA_ROOT_PATH UTOPIA_MODULES UTOPIA_HEADERS UTOPIA_SOURCES)

add_definitions(${UTOPIA_DEFS})
include_directories(${UTOPIA_INCLUDES})
include_directories(${UTOPIA_INCLUDE_PATHS})
include_directories(${UTOPIA_MODULES})
include_directories(${CMAKE_BINARY_DIR})
include_directories(.)

# link_libraries(${UTOPIA_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall ")# -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG ")

configure_file(${UTOPIA_ROOT_PATH}/utopia_Config.hpp.in ${CMAKE_BINARY_DIR}/utopia_Config.hpp)

add_library(utopia STATIC ${UTOPIA_SOURCES})
target_link_libraries(utopia ${UTOPIA_LIBRARIES})

if(WITH_CUDA)
	target_link_libraries(utopia utopia_cuda)
endif()

if(NOT UTOPIA_ARCHIVE_ONLY)
	add_library(utopia_dynamic SHARED ${UTOPIA_SOURCES})
	target_link_libraries(utopia_dynamic ${UTOPIA_LIBRARIES})
	
	if(WITH_CUDA)
		target_link_libraries(utopia_dynamic utopia_cuda)
	endif()

	add_executable(utopia_exec main.cpp ${UTOPIA_HEADERS})
	target_link_libraries(utopia_exec utopia)


	# target_link_libraries(utopia_exec utopia_dynamic)

	install(TARGETS utopia utopia_exec utopia_dynamic 
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
set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_INCLUDES=" "\n" "-I" "${UTOPIA_INCLUDES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_LIBRARIES=" "\n" "" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_makefile(CONFIG_CONTENT "UTOPIA_DEFS=" "\n" "" "${UTOPIA_DEFS}")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_FLAGS=${CMAKE_CXX_FLAGS}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_CXX_COMPILER=${CMAKE_CXX_COMPILER}\n")
set(CONFIG_CONTENT "${CONFIG_CONTENT}UTOPIA_C_COMPILER=${CMAKE_C_COMPILER}\n")


file(WRITE ${CMAKE_BINARY_DIR}/utopia_config.makefile ${CONFIG_CONTENT})

set(CONFIG_CONTENT "#Automatically generated by utopia CMakeLists.txt\n\n")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_INCLUDES\n" 		"\n)\n\n" "${UTOPIA_INCLUDES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_LIBRARIES\n" 	"\n)\n\n" "${UTOPIA_LIBRARIES}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_DEFS\n" 			"\n)\n\n" "${UTOPIA_DEFS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_FLAGS\n"  	"\n)\n\n" "${CMAKE_CXX_FLAGS}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_C_COMPILER\n"  	"\n)\n\n" "${CMAKE_C_COMPILER}")
append_list_to_string_for_cmake(CONFIG_CONTENT "LIST(APPEND UTOPIA_CXX_COMPILER\n" 	"\n)\n\n" "${CMAKE_CXX_COMPILER}")

file(WRITE ${CMAKE_BINARY_DIR}/utopia_config.cmake ${CONFIG_CONTENT})

install(FILES ${CMAKE_BINARY_DIR}/utopia_config.makefile DESTINATION config)
install(FILES ${CMAKE_BINARY_DIR}/utopia_config.cmake 	 DESTINATION config)

print_dependency_table()

