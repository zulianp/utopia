#include "utopia_SmoothApp.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/exodusII_io.h"

namespace utopia {

    void SmoothApp::run(Input &in)
    {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        in.get("mesh", mesh);

        MeshParamSmoother smoother;
        smoother.read(in);
        smoother.apply(mesh.mesh());

        libMesh::ExodusII_IO io(mesh.mesh());

        std::string output_path = "out.e";
        in.get("output-path", output_path);
        io.write(output_path);
    }
}
