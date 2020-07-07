set(UTOPIA_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

list(APPEND TEST_MODULES
    .
)

if(TARGET utopia_libmesh)
    list(APPEND TEST_MODULES libmesh)
endif()

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(UTOPIA_TEST_DIR TEST_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(utopia_fe_test PRIVATE ${LOCAL_SOURCES})
# utopia_link_default_targets(utopia_fe_test)

target_include_directories(utopia_fe_test PRIVATE ${UTOPIA_TEST_DIR})
target_include_directories(utopia_fe_test PRIVATE .)
target_include_directories(utopia_fe_test PRIVATE ${TEST_MODULES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_FE_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${UTOPIA_FE_DEV_FLAGS}")
