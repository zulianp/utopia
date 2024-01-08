# cmake_minimum_required(VERSION 2.8)

if(UTOPIA_ENABLE_MOONOLITH)

    set(UTOPIA_MOONOLITH_SEARCH_PATHS "${MOONOLITH_DIR};${MOONOLITH_DIR}/lib/cmake")

    if(UTOPIA_ENABLE_ENV_READ)
        set(UTOPIA_MOONOLITH_SEARCH_PATHS "${UTOPIA_MOONOLITH_SEARCH_PATHS};$ENV{MOONOLITH_DIR};$ENV{MOONOLITH_DIR}/lib/cmake")
    endif()

    if(UTOPIA_INSTALL_MOONOLITH)
        set(UTOPIA_MOONOLITH_SEARCH_PATHS "${MOONOLITH_INSTALL_DIR}")
        message(STATUS ${UTOPIA_MOONOLITH_SEARCH_PATHS})
    endif()

    find_package(
        ParMoonolith
        CONFIG
        HINTS
        ${UTOPIA_MOONOLITH_SEARCH_PATHS}
        PATH_SUFFIXES lib/cmake
        )

    if(ParMoonolith_FOUND)
        message(STATUS "Found ParMoonolith by config file")
        get_target_property(MOONOLITH_INCLUDES ParMoonolith::par_moonolith
                            INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(MOONOLITH_LIBRARIES ParMoonolith::par_moonolith
                            IMPORTED_LOCATION)
        if(NOT MOONOLITH_LIBRARIES)
            get_target_property(MOONOLITH_LIBRARIES ParMoonolith::par_moonolith
                                IMPORTED_LOCATION_RELEASE)
        endif()
        if(NOT MOONOLITH_LIBRARIES)
            get_target_property(MOONOLITH_LIBRARIES ParMoonolith::par_moonolith
                                IMPORTED_LOCATION_DEBUG)
        endif()
        return()
    else()
        message(WARNING "Could not find ParMoonolith by config file")
    endif()
endif()