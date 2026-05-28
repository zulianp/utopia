#ifndef UTOPIA_TRANSLATE_MESH_APP
#define UTOPIA_TRANSLATE_MESH_APP

#include "utopia_fe_Core.hpp"
#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"

namespace utopia {

    template <class Mesh>
    class TranslateMeshApp {
    public:
        static void run(Input &in) {
            std::string out;
            in.require("out", out);

            Mesh mesh;
            mesh.read(in);
            mesh.write(out);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_TRANSLATE_MESH_APP
