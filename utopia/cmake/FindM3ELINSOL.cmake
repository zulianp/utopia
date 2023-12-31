find_path(M3ELINSOL_INCLUDE_1
    NAMES "M3Elinsol_CXX.hpp"
    HINTS ../external/M3Elinsol
          ${M3ELINSOL_DIR}
          $ENV{M3ELINSOL_DIR}
    PATH_SUFFIXES CXX_ASPAMG_Driver
                  ../CXX_ASPAMG_Driver
    DOC "The M3ELINSOL_INCLUDES path"
)

find_path(M3ELINSOL_INCLUDE_2
    NAMES "M3Elinsol.hpp"
    HINTS ../external/M3Elinsol
          ${M3ELINSOL_DIR}
          ${M3ELINSOL_DIR}/../
          $ENV{M3ELINSOL_DIR}
          $ENV{M3ELINSOL_DIR}/../
    PATH_SUFFIXES Lib_src/INTERFACES
    DOC "The M3ELINSOL_INCLUDES path"
)

set(M3ELINSOL_INCLUDES "")
list(APPEND M3ELINSOL_INCLUDES
    ${M3ELINSOL_INCLUDE_1}
    ${M3ELINSOL_INCLUDE_2})



find_library(M3ELINSOL_LIBRARIES
    NAMES M3E_LINSOL
    PATHS ${M3ELINSOL_INCLUDES}/..
          ${M3ELINSOL_DIR}
          ${M3ELINSOL_DIR}/../
          $ENV{M3ELINSOL_DIR}
          $ENV{M3ELINSOL_DIR}/../
    PATH_SUFFIXES lib bin build
    DOC "The M3ELINSOL_LIBRARIES library to link against"
)


find_path(M3ELINSOL_SRC
    NAMES "M3Elinsol_CXX.cpp"
    HINTS ../external/M3Elinsol
          ${M3ELINSOL_INCLUDE_1}
          ${M3ELINSOL_INCLUDES}/..
          ${M3ELINSOL_DIR}
          $ENV{M3ELINSOL_DIR}
    PATH_SUFFIXES src
    DOC "The M3ELINSOL_SRC path"
)



list(APPEND
    M3ELINSOL_SOURCE_FILES
    "${M3ELINSOL_SRC}/M3Elinsol_CXX.cpp")

# message(STATUS "M3ELINSOL: ${M3ELINSOL_INCLUDES}, ${M3ELINSOL_LIBRARIES}, ${M3ELINSOL_SRC}")

message(STATUS "include: ${M3ELINSOL_INCLUDES}")
message(STATUS "lib: ${M3ELINSOL_LIBRARIES}")
message(STATUS "src: ${M3ELINSOL_SRC}")
message(STATUS "files: ${M3ELINSOL_SOURCE_FILES}")

if(M3ELINSOL_INCLUDES AND M3ELINSOL_LIBRARIES AND M3ELINSOL_SRC)
    set(M3ELINSOL_FOUND TRUE)
endif()
