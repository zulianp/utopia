
#include "utopia_DescribeFunctionSpaceApp.hpp"
#include "utopia_DescribeMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_stk.hpp"

void stk_describe_mesh(utopia::Input &in) {
    utopia::DescribeMeshApp<utopia::stk::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_describe_mesh);

void stk_describe_space(utopia::Input &in) {
    utopia::DescribeFunctionSpaceApp<utopia::stk::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_describe_space);
