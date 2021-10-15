#ifndef UTOPIA_CREATE_MESH_APP
#define UTOPIA_CREATE_MESH_APP

#include "utopia_fe_Core.hpp"
#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"

namespace utopia {

    template <class Mesh>
    class CreateMeshApp {
    public:
        static void run(Input &in) {
            Mesh mesh;
            mesh.read(in);

            std::string output = "out.e";
            in.get("output", output);
            mesh.write(output);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_CREATE_MESH_APP
