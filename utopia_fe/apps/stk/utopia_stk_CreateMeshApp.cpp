
#include "utopia_CreateMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_stk.hpp"

void stk_create_mesh(utopia::Input &in) {
    utopia::CreateMeshApp<utopia::stk::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_create_mesh);
