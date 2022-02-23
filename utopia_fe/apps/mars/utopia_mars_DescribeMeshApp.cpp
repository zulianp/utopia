
#include "utopia_DescribeMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_mars.hpp"

void mars_describe_mesh(utopia::Input &in) {
    utopia::DescribeMeshApp<utopia::mars::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(mars_describe_mesh);
