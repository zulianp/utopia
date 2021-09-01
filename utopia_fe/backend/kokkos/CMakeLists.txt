# TODO
utopiafe_add_module(utopia_fe_kokkos ${CMAKE_CURRENT_SOURCE_DIR} ".;fe;mesh;kernels;kernels/mech;materials;assembly;assembly/mech;assembly/phase_field;assembly/thermo_elasticity;utils;sampler;geom")
target_link_libraries(utopia_fe_kokkos PUBLIC utopia_fe_core Utopia::utopia)
