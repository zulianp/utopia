
#include "utopia_DescribeMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_stk_Mesh.hpp"

void stk_describe_mesh(utopia::Input &in) {
    utopia::DescribeMeshApp<utopia::stk::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_describe_mesh);
