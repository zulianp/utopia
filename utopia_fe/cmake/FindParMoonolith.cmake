cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

set(STAGE_DIR 				"${CMAKE_BINARY_DIR}/stage")
set(MOONOLITH_URL 			https://zulianp@bitbucket.org/zulianp/par_moonolith.git)
set(MOONOLITH_SOURCE_DIR 	${STAGE_DIR}/par_moonolith)
set(MOONOLITH_BIN_DIR 		${STAGE_DIR}/par_moonolith/bin)
set(MOONOLITH_INSTALL_DIR 	${CMAKE_BINARY_DIR}/external)

# set(MOONOLITH_INSTALL_DIR ${STAGE_DIR}/par_moonolith)

ExternalProject_Add(
	par_moonolith 
	PREFIX ${STAGE_DIR}
	GIT_REPOSITORY 		${MOONOLITH_URL}
	DOWNLOAD_DIR 		${STAGE_DIR}
	# SOURCE_DIR 		${MOONOLITH_SOURCE_DIR} 
	INSTALL_DIR         ${MOONOLITH_INSTALL_DIR}
	BINARY_DIR 			${MOONOLITH_SOURCE_DIR}
	# CONFIGURE_COMMAND   "${CMAKE_BINARY_DIR}/../cmake/install_par_moonolith.sh ${STAGE_DIR}/src/par_moonolith"
	# BUILD_COMMAND 		""
	# INSTALL_COMMAND 	""
	# BUILD_IN_SOURCE 	0   
	CMAKE_ARGS 			"-DCMAKE_INSTALL_PREFIX=${MOONOLITH_INSTALL_DIR}"
	LOG_CONFIGURE		1
	LOG_BUILD 			1
)

list(APPEND MOONOLITH_INCLUDES 
	${MOONOLITH_INSTALL_DIR}/include
	${MOONOLITH_INSTALL_DIR}/include/kernels
	)

set(MOONOLITH_LIB "")
list(APPEND MOONOLITH_LIB 
	"${MOONOLITH_INSTALL_DIR}/lib/libmoonolith_opencl.a"
	"${MOONOLITH_INSTALL_DIR}/lib/libpar_moonolith.a"
	"${MOONOLITH_INSTALL_DIR}/lib/libpar_moonolith_intersection.a"
	"${MOONOLITH_INSTALL_DIR}/lib/libpar_moonolith_mpi.a"
	"${MOONOLITH_INSTALL_DIR}/lib/libpar_moonolith_tree.a"
	"${MOONOLITH_INSTALL_DIR}/lib/libpar_moonolith_utils.a"
	)

