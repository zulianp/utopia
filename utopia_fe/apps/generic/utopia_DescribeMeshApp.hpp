#ifndef UTOPIA_DESCRIBE_MESH_APP
#define UTOPIA_DESCRIBE_MESH_APP

#include "utopia_fe_Core.hpp"
#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"

namespace utopia {

    template <class Mesh>
    class DescribeMeshApp {
    public:
        static void run(Input &in) {
            Mesh mesh;
            mesh.read(in);
            mesh.describe(utopia::out().stream());
        }
    };

}  // namespace utopia

#endif  // UTOPIA_DESCRIBE_MESH_APP
