find_package(
    xsimd
    CONFIG
    HINTS
    ${XSIMD_DIR}
    ${XSIMD_INCLUDES}
    $ENV{XSIMD_DIR}
    PATH_SUFFIXES
    lib/cmake)

if(xsimd_FOUND)
    target_include_directories(utopia SYSTEM PUBLIC ${xsimd_INCLUDE_DIRS})
    list(APPEND UTOPIA_THIRDPARTY_INCLUDES ${xsimd_INCLUDE_DIRS})
    set(UTOPIA_THIRDPARTY_INCLUDES
        ${UTOPIA_THIRDPARTY_INCLUDES}
        PARENT_SCOPE)

    set(UTOPIA_WITH_XSIMD TRUE)

else()
    set(UTOPIA_WITH_XSIMD FALSE)
endif()

set(UTOPIA_WITH_XSIMD
    {UTOPIA_WITH_XSIMD}
    PARENT_SCOPE)

list(APPEND XSIMD_MODULES .)
utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${XSIMD_MODULES}")
