#include "utopia_libmesh_Assembler.hpp"

// libmesh
#include "libmesh/mesh.h"

namespace utopia {

    bool libmesh_each(const LMFunctionSpace &space, std::function<bool(const libMesh::Elem &)> f) {
        auto &mesh = space.mesh().raw_type();

        auto e_begin = mesh.active_local_elements_begin();
        auto e_end = mesh.active_local_elements_end();

        bool ok = true;
        if (e_begin != e_end) {
            for (auto it = e_begin; it != e_end; ++it) {
                if (!f(**it)) {
                    ok = false;
                    break;
                }
            }
        }

        return ok;
    }

}  // namespace utopia