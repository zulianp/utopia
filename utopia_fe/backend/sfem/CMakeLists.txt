list(APPEND SFEM_MODULES .)

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${SFEM_MODULES}"
                   UTOPIA_FE_BUILD_INCLUDES UTOPIA_FE_HEADERS UTOPIA_FE_SOURCES)

  set(UTOPIA_FE_BUILD_INCLUDES
      ${UTOPIA_FE_BUILD_INCLUDES}
      PARENT_SCOPE)

  set(UTOPIA_FE_HEADERS
      ${UTOPIA_FE_HEADERS}
      PARENT_SCOPE)

  set(UTOPIA_FE_SOURCES
      ${UTOPIA_FE_SOURCES}
      PARENT_SCOPE)







find_path(
    SFEM_DIR
    NAMES  base/sfem_base.h
    HINTS ${SFEM_DIR} $ENV{SFEM_DIR}/sfem ${CMAKE_SOURCE_DIR}/external/sfem
)

if(SFEM_DIR)
    utopiafe_add_module(utopia_dep_sfem ${SFEM_DIR} 
        "algebra;base;baseline;eigen;matrix;mesh;operators;operators/cvfem/;operators/kernels/;operators/tet10/;operators/tet4/;operators/tri3/;operators/tri6/;operators/trishell3;operators/beam2/;perf;pizzastack;resampling;solver")

    target_link_libraries(utopia_dep_sfem PUBLIC Utopia::utopia)

    utopiafe_add_module(utopia_sfem ${CMAKE_CURRENT_SOURCE_DIR} ".")
    target_link_libraries(utopia_sfem PUBLIC utopia_dep_sfem utopia_fe_core)

    set(UTOPIA_FE_LIBRARIES
        ${UTOPIA_FE_LIBRARIES}
        PARENT_SCOPE)
else()
    message(FATAL_ERROR "Could not find SFEM!")
endif()
