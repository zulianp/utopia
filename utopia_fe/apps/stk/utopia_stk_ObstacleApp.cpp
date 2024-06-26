
#include "utopia_Main.hpp"

#include "utopia_ObstacleApp.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#ifdef UTOPIA_ENABLE_SFEM
#include "utopia_sfem_stk_SDFObstacle.hpp"
#endif  // UTOPIA_ENABLE_SFEM

void stk_obs(utopia::Input &in) {
#ifdef UTOPIA_ENABLE_SFEM
    utopia::register_sfem_stk_contact();
#endif  // UTOPIA_ENABLE_SFEM

    utopia::ObstacleApp<utopia::stk::FunctionSpace> obs;
    obs.run(in);
}

UTOPIA_REGISTER_APP(stk_obs);
