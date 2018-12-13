#First with try with clang compile or else
FIND_LIBRARY(MPI_CLANG_LIBRARY  
	NAMES mpi_cxx 
		  mpicxx
		  mpicxx-mpich-clang
	PATHS ${MPI_DIR}/lib
		  $ENV{MPI_DIR}/lib
		  /opt/local/lib/openmpi-mp/
		  /opt/local/lib/mpich-mp/
		  /opt/local/lib/mpich-clang/
		  /opt/local/lib

	DOC "The MPI_CLANG_LIBRARY library to link against"
) 


IF(MPI_CLANG_LIBRARY)
	SET(MPI_CXX_LIBRARIES ${MPI_CLANG_LIBRARY})

	GET_FILENAME_COMPONENT(MPI_LIB_DIR ${MPI_CLANG_LIBRARY} PATH)
	
	FIND_PATH(MPI_CLANG_HEADERS mpi.h
		HINTS ${MPI_DIR}/include
			  $ENV{MPI_DIR}/include
			  ${MPI_LIB_DIR}/../../include
			  ${MPI_LIB_DIR}/../include
			   /opt/local/include/openmpi-mp/
			  ${MPI_LIB_DIR}/../../include/openmpi-mp/
			  ${MPI_LIB_DIR}/../../include/mpich-clang
			  ${MPI_LIB_DIR}/../include/mpich-clang
			  /opt/local/include/mpich-clang
		DOC "The MPI_CLANG_HEADERS path"
	)	

	IF(MPI_CLANG_HEADERS)
		find_file(MPI_CXX_COMPILER 
			NAMES mpic++
			      mpicxx
			      mpicxx-openmpi-mp 
				  mpicxx-mpich-clang
			HINTS ${MPI_DIR}/bin
			      $ENV{MPI_DIR}/bin
			      ${MPI_CLANG_HEADERS}/../bin
			 	  ${MPI_CLANG_HEADERS}/../../bin
				  ${MPI_LIB_DIR}/../bin
				  ${MPI_LIB_DIR}/../../bin
				  /opt/local/bin/
			DOC "the MPI_CXX_COMPILER path"
		)

		find_file(MPI_C_COMPILER 
			NAMES mpicc-openmpi-mp
				  mpicc-mpich-clang
				  mpicc
			HINTS ${MPI_DIR}/bin
			      $ENV{MPI_DIR}/bin
			      ${MPI_CLANG_HEADERS}/../bin
			 	  ${MPI_CLANG_HEADERS}/../../bin
				  ${MPI_LIB_DIR}/../bin
				  ${MPI_LIB_DIR}/../../bin
				  /opt/local/bin/
			DOC "the MPI_C_COMPILER path"
		)

		IF(MPI_CXX_COMPILER AND MPI_C_COMPILER)
			#set variables
			SET(MPI_FOUND TRUE)
			SET(MPI_CXX_LIBRARIES ${MPI_CLANG_LIBRARY})
			SET(MPI_CXX_INCLUDE_PATH ${MPI_CLANG_HEADERS})
		ENDIF()
	ENDIF()
ENDIF()

# MESSAGE(STATUS "${MPI_CLANG_LIBRARY} ${MPI_CLANG_HEADERS} ${MPI_CXX_COMPILER}")


IF(NOT MPI_FOUND)
FIND_PACKAGE(MPI)
ENDIF()