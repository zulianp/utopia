# TODO
utopiafe_add_module(utopia_fe_kokkos ${CMAKE_CURRENT_SOURCE_DIR} ".;fe;mesh;kernels;kernels/mech;materials")
target_link_libraries(utopia_fe_kokkos PUBLIC utopia_fe_core Utopia::utopia)
