#include "utopia_Main.hpp"
#include "utopia_Multiphysics.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

void stk_multiphysics(utopia::Input &in) {
    utopia::Multiphysics<utopia::stk::FunctionSpace> multiphysics;
    multiphysics.read(in);
    multiphysics.run();
}

UTOPIA_REGISTER_APP(stk_multiphysics);
