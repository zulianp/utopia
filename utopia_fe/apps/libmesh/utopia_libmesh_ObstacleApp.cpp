
#include "utopia_Main.hpp"

#include "utopia_ObstacleApp.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"
#include "utopia_moonolith_libmesh_Obstacle.hpp"

void libmesh_obs(utopia::Input &in) {
    utopia::ObstacleProblem<utopia::libmesh::FunctionSpace> obs;
    obs.read(in);
    if (obs.valid) {
        obs.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(libmesh_obs);
