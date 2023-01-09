if(UTOPIA_ENABLE_VC)
  find_package(
    Vc
    CONFIG
    HINTS
    ${Vc_DIR}
    ${VC_DIR}
    ${VC_INCLUDES}
    $ENV{Vc_DIR}
    $ENV{VC_DIR}
    PATH_SUFFIXES
    lib/cmake)

  if(Vc_FOUND)
    # target_include_directories(utopia SYSTEM PUBLIC ${Vc_INCLUDE_DIR})
    # target_link_libraries(utopia PUBLIC ${Vc_LIBRARIES}) Why this does not
    # work? # target_link_libraries(utopia Vc::Vc)

    list(APPEND UTOPIA_BUILD_INCLUDES ${Vc_INCLUDE_DIRS})
    list(APPEND UTOPIA_DEP_LIBRARIES ${Vc_LIBRARIES})

    set(UTOPIA_BUILD_INCLUDES
        ${UTOPIA_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_DEP_LIBRARIES
        ${UTOPIA_DEP_LIBRARIES}
        PARENT_SCOPE)

    set(UTOPIA_WITH_VC TRUE)

    list(APPEND VC_MODULES .)
    # utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${VC_MODULES}")

    scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${VC_MODULES}"
                     UTOPIA_BUILD_INCLUDES UTOPIA_HEADERS UTOPIA_SOURCES)

    set(UTOPIA_BUILD_INCLUDES
        ${UTOPIA_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_HEADERS
        ${UTOPIA_HEADERS}
        PARENT_SCOPE)

    set(UTOPIA_SOURCES
        ${UTOPIA_SOURCES}
        PARENT_SCOPE)

    message(STATUS "${Vc_ARCHITECTURE_FLAGS}\n${Vc_ALL_FLAGS}")
  else()
    set(UTOPIA_WITH_VC FALSE)
  endif()
else()
  set(UTOPIA_WITH_VC FALSE)
endif()

set(UTOPIA_WITH_VC
    {UTOPIA_WITH_VC}
    PARENT_SCOPE)
