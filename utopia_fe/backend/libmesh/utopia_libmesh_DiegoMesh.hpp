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

        class Selector {
        public:
            std::set<int> blocks;
            std::vector<int> node_index;
            int n_elements;
            int n_nodes;

            void clear()
            {
                blocks.clear();
                node_index.clear();
                n_elements = 0;
                n_nodes = 0;
            }

        };


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

        void build_index(
            const libMesh::MeshBase &mesh,
            const std::set<int> &blocks,
            Selector &select)
        {
            select.blocks = blocks;
            build_index(mesh, blocks, select.n_elements, select.node_index, select.n_nodes);
        }

        bool write_coords(
            const Path &folder,
            const libMesh::MeshBase &mesh,
            const Selector &selector,
            const std::string &postfix = "");


        void build_index(
                const libMesh::MeshBase &mesh,
                const std::set<int> &blocks,
                int &n_elements,
                std::vector<int> &node_index,
                int &n_nodes);

        bool write_coords(
                const Path &folder,
                const libMesh::MeshBase &mesh,
                const std::vector<int> &node_index,
                const std::string &postfix = "");

        bool write(const Path &folder, const libMesh::MeshBase &mesh, const std::set<int> &blocks);
    };

}

#endif //UTOPIA_LIBMESH_DIEGOMESH_HPP
