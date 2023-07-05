#ifndef UTOPIA_RESCALE_MESH_APP
#define UTOPIA_RESCALE_MESH_APP

#include "utopia_fe_Core.hpp"
#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"

namespace utopia {

    template <class Mesh>
    class RescaleMeshApp {
    public:
        static void run(Input &in) {
            double scale = -1;
            in.require("scale", scale);

            std::string out;
            in.require("out", out);

            Mesh mesh;
            mesh.read(in);

            mesh.write(out);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_RESCALE_MESH_APP
