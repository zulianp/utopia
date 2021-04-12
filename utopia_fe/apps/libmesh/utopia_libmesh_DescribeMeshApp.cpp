
#include "utopia_DescribeMeshApp.hpp"
#include "utopia_Main.hpp"

#include "utopia_libmesh.hpp"

void libmesh_describe_mesh(utopia::Input &in) {
    utopia::DescribeMeshApp<utopia::libmesh::Mesh> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(libmesh_describe_mesh);
