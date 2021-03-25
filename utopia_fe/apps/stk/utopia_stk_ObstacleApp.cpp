
#include "utopia_Main.hpp"

#include "utopia_ObstacleApp.hpp"

#include "utopia_moonolith_stk_Obstacle.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

void stk_obs(utopia::Input &in) {
    utopia::ObstacleProblem<utopia::stk::FunctionSpace> obs;
    obs.read(in);
    if (obs.valid) {
        obs.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(stk_obs);
