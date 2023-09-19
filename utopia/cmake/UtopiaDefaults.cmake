# ##############################################################################
# ##############################################################################
# ##############################################################################
# FLAGS
if (WIN32)
set(UTOPIA_DEV_FLAGS
"-Wall"
)
else()

set(UTOPIA_DEV_FLAGS
    "-Wall -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder -Werror=return-type" # -Werror=uninitialized
)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=nonportable-include-path")
endif()

if(CMAKE_CXX_COMPILER MATCHES "nvcc")
  string(REPLACE "-Werror=reorder" "" UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS}")
  message(
    STATUS
      "Using nvcc so we remove the -Werror-reorder flag due to a bug in nvcc")
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wno-clobbered -Wno-vla -Wno-pragmas
# -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-literal-suffix
# -Wno-deprecated-declarations -Wno-misleading-indentation
# -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-nonnull-compare
# -Wno-address -Wno-inline -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS -Werror"
# CACHE STRING "Warnings as errors settings")

# set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS}
# -Werror=inconsistent-missing-override")

if(UTOPIA_ENABLE_CLANG_TIDY)
  # FIXME
  find_program(
    CLANG_TIDY
    NAMES clang-tidy-mp-7.0 clang-tidy-7 clang-tidy-6.0 clang-tidy-5.0
          clang-tidy-4.0 clang-tidy
    PATHS /opt/local/bin)

  set(CMAKE_CXX_CLANG_TIDY
      "${CLANG_TIDY};-format-style=file;-fix;-export-fixes=fixes.yml"
  )# -warnings-as-errors=*; #-header-filter=.
endif()

# More annoying
if(NOT WIN32)
set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Wextra ")
endif()
# More restrictive
if(UTOPIA_PULL_REQUEST_MODE)
  set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=deprecated-declarations ")
  # set(UTOPIA_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Werror=unused-variable
  # -Werror=unused-local-typedef ")
endif()

# include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/external/GSL/include)



if(WIN32)
  set(CMAKE_CXX_FLAGS_DEBUG
      "${CMAKE_CXX_FLAGS_DEBUG}   -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
endif()

message(
  STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if(UTOPIA_ENABLE_SANITIZER)
  # if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")

  set(UTOPIA_DEV_FLAGS
      "${UTOPIA_DEV_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls "
  )

  # endif()
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################