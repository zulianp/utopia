
#include "utopia_Main.hpp"

#include "utopia_ObstacleApp.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"
#include "utopia_moonolith_libmesh_Obstacle.hpp"

void libmesh_obs(utopia::Input &in) {
    utopia::ObstacleApp<utopia::libmesh::FunctionSpace> obs;
    obs.run(in)
}

UTOPIA_REGISTER_APP(libmesh_obs);
