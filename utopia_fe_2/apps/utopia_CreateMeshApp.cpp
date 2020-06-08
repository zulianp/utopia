#include "utopia_AppRunner.hpp"
#include "utopia_libmesh.hpp"

#include <cmath>

namespace utopia {

    static void create_mesh(Input &in) {
        using Comm = Traits<LMMesh>::Communicator;

        Comm comm;
        LMMesh mesh(comm);
        mesh.read(in);

        std::string path = "out.e";
        in.get("output", path);

        mesh.write(path);
    }

    UTOPIA_REGISTER_APP(create_mesh);

}  // namespace utopia