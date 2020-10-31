cmake_minimum_required(VERSION 2.8)

if(MOONOLITH_DIR OR DEFINED ENV{MOONOLITH_DIR})

    find_package(
        ParMoonolith
        CONFIG
        HINTS
        ${MOONOLITH_DIR}
        ${MOONOLITH_DIR}/lib/cmake
        $ENV{MOONOLITH_DIR}
        $ENV{MOONOLITH_DIR}/lib/cmake)

    if(ParMoonolith_FOUND)
        message(STATUS "Found ParMoonolith by config file")
        return()
    else()
        message(FATAL_ERROR "Could not find ParMoonolith by config file")
    endif()

    find_path(MOONOLITH_INSTALLATION_PATH NAME config/moonolith_config.cmake
              HINTS ${MOONOLITH_DIR} $ENV{MOONOLITH_DIR})

    if(MOONOLITH_INSTALLATION_PATH)
        message(
            STATUS
                "Found moonolith installation at ${MOONOLITH_INSTALLATION_PATH}"
        )
        include(${MOONOLITH_INSTALLATION_PATH}/config/moonolith_config.cmake)
        include(FindPackageHandleStandardArgs)

        find_package_handle_standard_args(
            MOONOLITH REQUIRED_VARS MOONOLITH_LIBRARIES MOONOLITH_INCLUDES)

        mark_as_advanced(MOONOLITH_INCLUDES MOONOLITH_LIBRARIES)

        if(MOONOLITH_FOUND)
            add_custom_target(par_moonolith)
        endif()
    else()

    endif()

endif()

if(NOT MOONOLITH_FOUND OR FORCE_INSTALL_MOONOLITH)
    include(FetchContent)
    message(STATUS "Fetching par_moonolith, since it could not be found.")

    set(MOONOLITH_ENABLE_BENCHMARK
        OFF}
        CACHE INTERNAL "")

    set(MOONOLITH_ENABLE_TESTING
        OFF
        CACHE INTERNAL "")

    FetchContent_Declare(
        moonolith
        GIT_REPOSITORY https://bitbucket.org/zulianp/par_moonolith.git)
    FetchContent_MakeAvailable(moonolith)

    add_library(ParMoonolith::par_moonolith ALIAS par_moonolith)
endif()
