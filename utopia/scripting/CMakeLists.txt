list(APPEND SCRIPTING_MODULES
    .
)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR SCRIPTING_MODULES LOCAL_HEADERS LOCAL_SOURCES)

add_library(utopia_c SHARED ${LOCAL_SOURCES})
target_link_libraries(utopia_c utopia)
target_include_directories(utopia_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


install(TARGETS utopia_c
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


find_package(SWIG REQUIRED)

INCLUDE(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_SWIG_FLAGS "")
set_property(SOURCE utopia.i PROPERTY CPLUSPLUS ON)
swig_add_library(utopia TYPE SHARED LANGUAGE python SOURCES utopia.i)
swig_link_libraries(utopia ${PYTHON_LIBRARIES} utopia_c)

