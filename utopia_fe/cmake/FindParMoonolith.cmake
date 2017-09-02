cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

set(STAGE_DIR 				"${CMAKE_BINARY_DIR}/stage")
set(MOONOLITH_URL 			https://zulianp@bitbucket.org/zulianp/par_moonolith.git)
set(MOONOLITH_SOURCE_DIR 	${STAGE_DIR}/par_moonolith)
set(MOONOLITH_BIN_DIR 		${STAGE_DIR}/par_moonolith/bin)
set(MOONOLITH_INSTALL_DIR 	${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(
	par_moonolith 
	UPDATE_COMMAND		"" #FIXME
	PREFIX ${STAGE_DIR}
	GIT_REPOSITORY 		${MOONOLITH_URL}
	DOWNLOAD_DIR 		${STAGE_DIR} 
	INSTALL_DIR         ${MOONOLITH_INSTALL_DIR}
	BINARY_DIR 			${MOONOLITH_SOURCE_DIR}
	CMAKE_ARGS 			"-DCMAKE_INSTALL_PREFIX=${MOONOLITH_INSTALL_DIR};-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE};-DENABLE_SANITIZER=${ENABLE_SANITIZER}"
	LOG_CONFIGURE		1
	LOG_BUILD 			1
)

list(APPEND MOONOLITH_INCLUDES 
	${MOONOLITH_INSTALL_DIR}/include
	${MOONOLITH_INSTALL_DIR}/include/kernels
	)

set(MOONOLITH_LIB "")
list(APPEND MOONOLITH_LIB 
	"-L${MOONOLITH_INSTALL_DIR}/lib"
	"-lmoonolith_opencl"
	"-lpar_moonolith"
	"-lpar_moonolith_intersection"
	"-lpar_moonolith_mpi"
	"-lpar_moonolith_tree"
	"-lpar_moonolith_utils"
	)
