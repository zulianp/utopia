
#include "utopia_Main.hpp"
#include "utopia_libmesh.hpp"

using namespace utopia;

void libmesh_refine(utopia::Input &in) {
    int n = 1;
    Path output = "./out.e";
    utopia::libmesh::Mesh mesh;
    mesh.read(in);

    if (mesh.empty()) {
        utopia::err() << "[Error] failed to initialize mesh\n";
        return;
    }

    in.get("output", output);
    in.get("n", n);

    mesh.uniform_refine(n);

    if (!mesh.write(output)) {
        utopia::err() << "[Error] failed to write mesh at " << output << "\n";
    }
}

UTOPIA_REGISTER_APP(libmesh_refine);
