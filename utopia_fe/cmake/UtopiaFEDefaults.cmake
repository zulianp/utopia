set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
          FORCE)

set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING
          "Flags used by the C++ compiler during AddressSanitizer builds."
          FORCE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

# To kill the policy warning  (maybe not a good idea yet)
set(CMAKE_MACOSX_RPATH 1)


if(UTOPIA_ENABLE_TRILINOS OR UTOPIA_ENABLE_MARS OR UTOPIA_ENABLE_INTREPID2)
    set(UTOPIA_ENABLE_KOKKOS TRUE)
endif()


if(NOT UTOPIA_ENABLE_MARS)
    if(UTOPIA_ENABLE_WARNINGS)

        set(UTOPIA_FE_DEV_FLAGS
            "-Wall -Werror=enum-compare -Werror=delete-non-virtual-dtor -Werror=reorder -Werror=return-type" # -Werror=uninitialized
)

        set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_DEV_FLAGS} -Wextra ")

    else()
        set(UTOPIA_FE_DEV_FLAGS "-w -Wno-deprecated-declarations")
    endif()

endif()

# More annoying set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_FE_DEV_FLAGS} -Wextra
# -Werror=uninitialized ")

# More restrictive set(UTOPIA_FE_DEV_FLAGS "${UTOPIA_FE_DEV_FLAGS}
# -Werror=unused-variable -Werror=unused-local-typedef
# -Werror=deprecated-declarations ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(UTOPIA_ENABLE_TRILINOS_ALGEBRA)
    set(UTOPIA_ENABLE_TRILINOS_ALGEBRA TRUE)
endif()

if(UTOPIA_ENABLE_NEW_TRANSFER)
    set(UTOPIA_ENABLE_NEW_TRANSFER TRUE)
endif()

set(UTOPIA_FE_ROOT_PATH ${CMAKE_SOURCE_DIR})
set(UTOPIA_FE_CMAKES_PATH ${UTOPIA_FE_ROOT_PATH}/libcmake)

if(UTOPIA_ENABLE_SANITIZER AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1")
        set(CMAKE_CXX_FLAGS_DEBUG
            "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls "
        )
    endif()
endif()