if(Trilinos_Kokkos_FOUND)

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
