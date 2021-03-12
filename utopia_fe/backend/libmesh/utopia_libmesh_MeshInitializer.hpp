#ifndef UTOPIA_LIBMESH_MESH_INITIALIZER_HPP
#define UTOPIA_LIBMESH_MESH_INITIALIZER_HPP

#include "utopia_Input.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"

namespace utopia {
    namespace libmesh {

        class MeshInitializer final : public Configurable {
        public:
            MeshInitializer(Mesh &mesh);
            void read(Input &in) override;

        private:
            Mesh &mesh_;

            class Build;
            class Shift;
            class Scale;
        };

    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_MESH_INITIALIZER_HPP
