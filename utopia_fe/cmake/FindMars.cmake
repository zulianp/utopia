if(UTOPIA_ENABLE_MARS)
  set(MARS_SEARCH_PATHS "/usr/local/;/usr/;${MARS_DIR}
        ${MARS_INCLUDES}")

  if(UTOPIA_ENABLE_ENV_READ)
    set(MARS_SEARCH_PATHS "${MARS_SEARCH_PATHS};$ENV{MARS_DIR};$ENV{Mars_DIR}")
  endif()

  find_package(Mars CONFIG HINTS ${MARS_SEARCH_PATHS} PATH_SUFFIXES lib/cmake)

  if(Mars_FOUND)
    message(STATUS "Found Mars by config file")

    get_target_property(MARS_INCLUDES Mars::mars INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(MARS_LIBRARIES Mars::mars IMPORTED_LOCATION)

    if(NOT MARS_LIBRARIES)
      get_target_property(MARS_LIBRARIES Mars::mars IMPORTED_LOCATION_RELEASE)
    endif()

    if(NOT MARS_LIBRARIES)
      get_target_property(MARS_LIBRARIES Mars::mars IMPORTED_LOCATION_DEBUG)
    endif()

    if(MARS_ENABLE_ADIOS)
      get_target_property(ADIOS_LIBRARIES adios2::adios2
                          INTERFACE_LINK_LIBRARIES)
    endif()

    # if(MARS_ENABLE_VTK) 
    #   get_target_property(VTK_LIBRARIES adios2:: INTERFACE_LINK_LIBRARIES)
    # endif()
  endif()

endif()
