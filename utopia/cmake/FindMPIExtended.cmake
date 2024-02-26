# First with try with clang compile

if(TPL_MPI_LIBRARIES AND TPL_MPI_INCLUDE_DIRS)
  set(MPI_FOUND TRUE)
  set(MPI_CXX_LIBRARIES ${TPL_MPI_LIBRARIES})
  set(MPI_CXX_INCLUDE_PATH ${TPL_MPI_INCLUDE_DIRS})
  return()
endif()

if(MPI_CXX_INCLUDE_PATH AND MPI_CXX_LIBRARIES)
  set(MPI_FOUND TRUE)
  return()
endif()


if(UTOPIA_ENABLE_ENV_READ)
if(NOT MPI_DIR)
  set(MPI_DIR $ENV{MPI_DIR})
endif()
endif()

if(${MPI_DIR})
  set(MPI_SEARCH_PATHS_LIBRARY
      "${MPI_DIR}/lib;${MPI_LIB_DIR};/opt/local/lib/openmpi-mp/;/opt/local/lib/mpich-mp/;/opt/local/lib/mpich-clang/;/opt/local/lib"
  )
endif()

if(${MPI_DIR})
  set(MPI_SEARCH_PATHS_HEADERS
      "${MPI_DIR}/include;${MPI_INCLUDE_DIR};${MPI_TEMP_LIBRARY}/../include;${MPI_LIB_DIR}/../../include;${MPI_LIB_DIR}/../include;/opt/local/include/openmpi-mp/;${MPI_LIB_DIR}/../../include/openmpi-mp/;${MPI_LIB_DIR}/../../include/mpich-clang;${MPI_LIB_DIR}/../include/mpich-clang;/opt/local/include/mpich-clang"
  )
endif()

if(UTOPIA_ENABLE_ENV_READ)
  set(MPI_SEARCH_PATHS_LIBRARY "${MPI_SEARCH_PATHS_LIBRARY};$ENV{MPI_DIR}/lib")
  set(MPI_SEARCH_PATHS_HEADERS
      "${MPI_SEARCH_PATHS_HEADERS};$ENV{MPI_DIR}/include;$ENV{MPI_INCLUDE_DIR}")
endif()

if(APPLE)

  find_library(
    MPI_TEMP_LIBRARY
    NAMES mpi_cxx mpicxx-mpich-clang
    PATHS ${MPI_SEARCH_PATHS_LIBRARY}
    DOC "The MPI_TEMP_LIBRARY library to link against")

  if(MPI_TEMP_LIBRARY)
    set(MPI_CXX_LIBRARIES ${MPI_TEMP_LIBRARY})

    get_filename_component(MPI_LIB_DIR ${MPI_TEMP_LIBRARY} PATH)

    find_path(
      MPI_TEMP_HEADERS mpi.h
      HINTS ${MPI_SEARCH_PATHS_HEADERS}
      DOC "The MPI_TEMP_HEADERS path")

    if(MPI_TEMP_HEADERS)
      find_file(
        MPI_CXX_COMPILER
        NAMES mpic++ mpicxx-openmpi-mp mpicxx-mpich-clang
        HINTS ${MPI_TEMP_HEADERS}/../bin ${MPI_TEMP_HEADERS}/../../bin
              ${MPI_LIB_DIR}/../bin ${MPI_LIB_DIR}/../../bin /opt/local/bin/
        DOC "the MPI_COMPILER_PATH dir")

      find_file(
        MPI_C_COMPILER
        NAMES mpicc mpicc-openmpi-mp mpicc-mpich-clang
        HINTS ${MPI_TEMP_HEADERS}/../bin ${MPI_TEMP_HEADERS}/../../bin
              ${MPI_LIB_DIR}/../bin ${MPI_LIB_DIR}/../../bin /opt/local/bin/
        DOC "the MPI_COMPILER_PATH dir")

      if(MPI_CXX_COMPILER AND MPI_C_COMPILER)
        # set variables
        set(MPI_FOUND TRUE)
        set(MPI_CXX_LIBRARIES ${MPI_TEMP_LIBRARY})
        set(MPI_CXX_INCLUDE_PATH ${MPI_TEMP_HEADERS})
      endif()
    endif()
  endif()
endif(APPLE)

# MESSAGE(STATUS "${MPI_TEMP_LIBRARY} ${MPI_TEMP_HEADERS} ${MPI_CXX_COMPILER}")

if(NOT MPI_FOUND)
  find_package(MPI)
  set(MPIExtended_FOUND ${MPI_FOUND})
endif()
