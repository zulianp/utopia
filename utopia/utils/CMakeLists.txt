list(APPEND UTILS_MODULES
    .
    action
    unit_testing
    app_management
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR UTILS_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_utils STATIC ${LOCAL_SOURCES})
utopia_link_default_targets(utopia_utils)

target_include_directories(utopia_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})
target_include_directories(utopia_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(utopia_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/unit_testing)
target_include_directories(utopia_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/action)
target_include_directories(utopia_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/app_management)

find_package(MPIExtended)

if(TRY_WITH_MPI)
    if(MPI_FOUND)
        list(APPEND UTOPIA_INCLUDES ${MPI_CXX_INCLUDE_PATH})
        list(APPEND UTOPIA_LIBRARIES  ${MPI_CXX_LIBRARIES} )

        target_include_directories(utopia_utils SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
        # target_link_libraries(utopia_utils utopia_core ${MPI_CXX_LIBRARIES} )
        target_link_libraries(utopia_utils ${MPI_CXX_LIBRARIES})


        message(STATUS "WHAT: ${MPI_CXX_INCLUDE_PATH} ${MPI_CXX_LIBRARIES}")

        set(WITH_MPI TRUE PARENT_SCOPE)
        if(MPI_DIR)
            set(MPI_DIR ${MPI_DIR} PARENT_SCOPE)
        endif()
    endif()
else()
    message(WARNING "NO Proper MPI installation")
endif()



set_utopia_compiler_features(utopia_utils)

install(TARGETS utopia_utils
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_utils" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")


set(UTOPIA_INCLUDES ${UTOPIA_INCLUDES} PARENT_SCOPE)
set(UTOPIA_LIBRARIES ${UTOPIA_LIBRARIES} PARENT_SCOPE)
set(UTOPIA_DEFS ${UTOPIA_DEFS} PARENT_SCOPE)
