# TODO
utopiafe_add_module(utopia_fe_kokkos ${CMAKE_CURRENT_SOURCE_DIR} ".;fe;fe/lagrange;mesh;kernels;kernels/mech;kernels/thermo_elasticity;materials;assembly;assembly/generated;assembly/mech;assembly/mech/generated;assembly;assembly/mech;assembly/mech/generated/2D;assembly;assembly/mech;assembly/mech/generated/3D;assembly/phase_field;assembly/thermo_elasticity;utils;sampler;geom;fe/generated")
target_link_libraries(utopia_fe_kokkos PUBLIC utopia_fe_core Utopia::utopia)


set(UTOPIA_FE_LIBRARIES
    ${UTOPIA_FE_LIBRARIES}
    PARENT_SCOPE)
