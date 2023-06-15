

#include "utopia_Main.hpp"
#include "utopia_RescaleMeshApp.hpp"

#include "utopia_stk.hpp"

void stk_rescale_mesh(utopia::Input &in) {
    utopia::RescaleMeshApp<utopia::stk::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_rescale_mesh);
