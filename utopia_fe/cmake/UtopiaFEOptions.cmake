option(UTOPIA_ENABLE_SANITIZER "check for memory access problems" OFF)
option(UTOPIA_ENABLE_BOOST "enable boost features" OFF)
option(UTOPIA_ENABLE_MOONOLITH_PROFILER
       "enable moonolith profiling capabilities" OFF)

option(UTOPIA_ENABLE_CXX14_FEATURES "Enable usage of cxx14 standard" ON)
option(UTOPIA_ENABLE_CXX17_FEATURES "Enable usage of cxx17 standard" ON)

option(UTOPIA_ENABLE_MOOSE_ENV_MODE
       "Allows to use the moose installation compilers" ON)
option(UTOPIA_ENABLE_TRILINOS_ALGEBRA
       "Allows to use the moose installation compilers" OFF)
option(UTOPIA_ENABLE_MARS "Tries to use mars backend" OFF)
# option(UTOPIA_ENABLE_MARS_VTK "Enable vtk output for mars" OFF)
# option(UTOPIA_ENABLE_MARS_ADIOS2 "Enable adios2 output for mars" OFF)
option(UTOPIA_ENABLE_INTREPID2 "Tries to use the intrepid2 related features" OFF)
option(UTOPIA_ENABLE_NEW_TRANSFER "Use new transfer features" ON)
# option(UTOPIA_INSTALL_LIBMESH "install libmesh" OFF)
option(UTOPIA_INSTALL_MOONOLITH "Local install of par_moonolith." OFF)
option(UTOPIA_ENABLE_WARNINGS "Compiler warnings" ON)
option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)

option(UTOPIA_ENABLE_LIBMESH "Enable the libmesh functionalities" OFF)
option(UTOPIA_ENABLE_LIBMESH_LEGACY "Enable the legacy utopia::libmesh functionalities" OFF)
option(UTOPIA_ENABLE_LIBMESH_KOKKOS "Enable the libmesh functionalities" OFF)

option(UTOPIA_INSTALL_LIBMESH "Local install of libmesh." OFF)

# option(UTOPIA_ENABLE_PETSC "Enable the PETSc DM functionalities" OFF)
option(UTOPIA_ENABLE_PETSCDM "Enable the PETSc DM functionalities" OFF)

# Once we have cleaner separation between this and legacy turn OFF
option(UTOPIA_ENABLE_LIBMESH_DEPRECATED "Enable the deprecated utopia::libmesh functionalities" OFF)


option(UTOPIA_ENABLE_TESTS "Enable utopia_fe testing." ON)

option(UTOPIA_ENABLE_STK "Enable the stk functionalities" OFF)
option(UTOPIA_ENABLE_MOONOLITH "Enable the moonolith functionalities" ON)
option(UTOPIA_ENABLE_ARBORX "Enable the ArborX backend" OFF)

option(UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE "Removed trilinos deprecated code" OFF)

option(UTOPIA_ENABLE_FLUYA_MODE "Create utopia_fe configuration required by Fluya (overrides everything else)" OFF)
option(UTOPIA_ENABLE_SFEM "Enable SFEM backend" OFF)


option(UTOPIA_ENABLE_ENV_READ "Enable the installation to look for env variables." ON)


set(UTOPIA_FE_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTOPIA_FE_ROOT_PATH}/cmake")