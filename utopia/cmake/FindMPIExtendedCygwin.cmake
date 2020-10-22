# Check if Cygwin is installed (windows system) if so do not use CLANG
    message(STATUS "On windows, skipping preamble for FindMPIExtended")

    find_library(MPI_C_LIBRARY
    NAMES mpi
    PATHS ${MPI_DIR}/lib
          $ENV{MPI_DIR}/lib
          /usr/lib
          DOC "The MPI_C_LIBRARY library to link against"
    )

    if(NOT MPI_C_LIBRARY)
      message(FATAL_ERROR "Could not find MPI_C_LIBRARY")
    endif()

    list(APPEND MPI_CXX_LIBRARIES ${MPI_C_LIBRARY})

    get_filename_component(MPI_LIB_DIR ${MPI_C_LIBRARY} PATH)

        find_path(MPI_C_HEADERS mpi.h
        HINTS ${MPI_DIR}/include
              $ENV{MPI_DIR}/include
              ${MPI_LIB_DIR}/../../include
              ${MPI_LIB_DIR}/../include
              /usr/include
        DOC "The MPI_C_HEADERS path"
    )

    #Look for excecutables of CXX_Compiler
    IF(MPI_C_HEADERS)
        find_file(MPI_CXX_COMPILER
            NAMES mpic++
                  mpicxx
                  mpicxx-openmpi-mp
                  mpicxx-mpich-clang
                  mpic++-mpich-clang
                  mpic++-mpich-mp
            HINTS ${MPI_DIR}/bin
                  $ENV{MPI_DIR}/bin
                  ${MPI_C_HEADERS}/../bin
                  /usr/bin
            DOC "the MPI_CXX_COMPILER path"
        )

        #Look for executables of C_Compiler
        find_file(MPI_C_COMPILER
            NAMES mpicc-openmpi-mp
                  mpicc-mpich-clang
                  mpicc-mpich-mp
                  mpicc
            HINTS ${MPI_DIR}/bin
                  $ENV{MPI_DIR}/bin
                  ${MPI_C_HEADERS}/../bin
                  /usr/bin
            DOC "the MPI_C_COMPILER path"
        )

          set(MPI_FOUND TRUE)

            list(APPEND MPI_CXX_INCLUDE_PATH ${MPI_C_HEADERS})
            message(STATUS "MPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}\n-- MPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}\n-- MPI_CXX_COMPILER=${MPI_CXX_COMPILER}\n-- MPI_C_COMPILER=${MPI_C_COMPILER}")
    endif()