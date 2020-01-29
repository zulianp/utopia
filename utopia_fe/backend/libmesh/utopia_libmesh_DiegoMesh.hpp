#ifndef UTOPIA_LIBMESH_DIEGOMESH_HPP
#define UTOPIA_LIBMESH_DIEGOMESH_HPP

#include "utopia_Path.hpp"
#include "utopia_Input.hpp"
#include "utopia_fe_base.hpp"

#include <string>

namespace libMesh {
    class MeshBase;
}

namespace utopia {

    class DiegoMeshWriter : public Configurable {
    public:
        void read(Input &in) override;
        bool write_field(
            const Path &folder,
            const std::string &field_name,
            const libMesh::MeshBase &mesh,
            const UVector &field,
            const int component,
            const int sys_num = 0);
        bool write_headers(const Path &folder, const libMesh::MeshBase &mesh);
        bool write_elements(const Path &folder, const libMesh::MeshBase &mesh);
        bool write_coords(const Path &folder, const libMesh::MeshBase &mesh, const std::string &postfix = "");
        bool write(const Path &folder, const libMesh::MeshBase &mesh);
    };

}

#endif //UTOPIA_LIBMESH_DIEGOMESH_HPP
