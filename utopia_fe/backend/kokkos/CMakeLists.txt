if(Kokkos_FOUND)

  list(
    APPEND
    KOKKOS_MODULES
    .
    assembly
    assembly/flow
    assembly/generated
    assembly/mech
    assembly/mech/old
    assembly/mech/generated
    assembly/mech/generated/2D
    assembly/mech/generated/3D
    assembly/old
    assembly/phase_field
    assembly/phase_field/old
    assembly/thermo_elasticity
    assembly/thermo_elasticity/old
    fe
    fe/generated
    fe/lagrange

    geom
    kernels
    kernels/mech
    kernels/thermo_elasticity
    mesh
    sampler
    utils)
  message(STATUS "KOKKOS_MODULES:${KOKKOS_MODULES}")

  scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${KOKKOS_MODULES}"
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

endif()

# find_package( Kokkos NO_MODULE NO_CMAKE_ENVIRONMENT_PATH
# NO_SYSTEM_ENVIRONMENT_PATH HINTS $ENV{KOKKOS_DIR}
# $ENV{TRILINOS_DIR}/lib/cmake/Kokkos ${Trilinos_DIR}/../Kokkos REQUIRED)

# utopiafe_add_module( utopia_fe_kokkos ${CMAKE_CURRENT_SOURCE_DIR}
# ".;fe;fe/lagrange;mesh;kernels;kernels/mech;kernels/mech/old;kernels/thermo_elasticity;kernels/thermo_elasticity/old;materials;assembly;assembly/generated;assembly/mech;assembly/mech/generated;assembly;assembly/old;assembly/mech;assembly/mech/old;assembly/mech/generated/2D;assembly;assembly/mech;assembly/mech/generated/3D;assembly/phase_field;assembly/phase_field/old;assembly/thermo_elasticity;assembly/thermo_elasticity/old;utils;sampler;geom;fe/generated;assembly/flow"
# ) target_link_libraries(utopia_fe_kokkos PUBLIC utopia_fe_core Utopia::utopia)
# target_link_libraries(utopia_fe_kokkos PUBLIC Kokkos::kokkos)
# target_include_directories(utopia_fe_kokkos PUBLIC ${Kokkos_INCLUDE_DIRS})

# set(UTOPIA_FE_LIBRARIES ${UTOPIA_FE_LIBRARIES} PARENT_SCOPE)
