# apps.cmake

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
set(APPS_MODULES "")

set(UTOPIA_APPS_DIR ${CMAKE_SOURCE_DIR}/apps)

list(APPEND APPS_MODULES generic)

if(UTOPIA_ENABLE_LIBMESH)
    list(APPEND APPS_MODULES libmesh utopia_fe)
endif()

if(UTOPIA_ENABLE_STK AND UTOPIA_ENABLE_INTREPID2)
    list(APPEND APPS_MODULES stk intrepid2)
endif()

find_project_files(${UTOPIA_APPS_DIR} "${APPS_MODULES}" LOCAL_HEADERS
                   LOCAL_SOURCES)
target_sources(
    utopia_fe_exec
    PRIVATE ${LOCAL_SOURCES}
    PRIVATE ${LOCAL_HEADERS})

target_include_directories(utopia_fe_exec PRIVATE ${UTOPIA_APPS_DIR}
                                                  ${UTOPIA_APPS_DIR}/generic)

target_include_directories(utopia_fe_exec PRIVATE .)

foreach(MODULE ${APPS_MODULES})
    target_include_directories(utopia_fe_exec
                               PRIVATE ${UTOPIA_APPS_DIR}/${MODULE})
endforeach()
