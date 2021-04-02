#include "utopia_Main.hpp"
#include "utopia_Multiphysics.hpp"

#include "utopia_libmesh_OmniAssembler.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"
#include "utopia_moonolith_libmesh_Obstacle.hpp"

void libmesh_multiphysics(utopia::Input &in) {
    utopia::Multiphysics<utopia::libmesh::FunctionSpace> multiphysics;
    multiphysics.read(in);
    multiphysics.run();
}

UTOPIA_REGISTER_APP(libmesh_multiphysics);
