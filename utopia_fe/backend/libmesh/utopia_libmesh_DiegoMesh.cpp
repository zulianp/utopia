#include "utopia_libmesh_DiegoMesh.hpp"

#include <fstream>

#include "libmesh/mesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/utility.h"
#include "libmesh/enum_elem_type.h"



namespace utopia {

    void DiegoMeshWriter::read(Input &in)
    {

    }

    bool DiegoMeshWriter::write_field(
        const Path &folder,
        const std::string &field_name,
        const libMesh::MeshBase &mesh,
        const UVector &field,
        const int component,
        const int sys_num)
    {

        std::string path = folder / (field_name + ".raw");

        std::ofstream os;

        os.open(path);
        if(!os.good()) {
            std::cout << "no file at " << path << std::endl;
            return false;
        }

        Read<UVector> r_d(field);

        auto m_it  = mesh.local_nodes_begin();
        auto m_end = mesh.local_nodes_end();

        for(; m_it != m_end; ++m_it) {
            const int dof_id = (*m_it)->dof_number(sys_num, component, 0);
            float v = field.get(dof_id);
            os.write((char *)&v, sizeof(v));
        }

        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_headers(const Path &folder, const libMesh::MeshBase &mesh)
    {
        std::string node_count_path = folder / "node_count.txt";
        std::string element_count_path = folder / "element_count.txt";
        std::string format_file = folder / "format.txt";

        std::ofstream os;

        // os.open(format_file);
        // if(!os.good()) {
        //     std::cout << "no file at " << format_file << std::endl;
        //     return false;
        // }

        // os << libMesh::Utility::enum_to_string((*mesh.active_local_elements_begin())->type()) << "\n";

        // os.close();

        os.open(node_count_path.c_str());

        if(!os.good()) {
            std::cout << "no file at " << node_count_path << std::endl;
            return false;
        }

        os << mesh.n_nodes() << "\n";
        os.close();

        os.open(element_count_path.c_str());
        if(!os.good()) {
            std::cout << "no file at " << element_count_path << std::endl;
            return false;
        }

        os << mesh.n_active_elem() << "\n";
        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_elements(const Path &folder, const libMesh::MeshBase &mesh)
    {
        std::string elem_file  = folder / "elem.raw";
        std::ofstream os;
        os.open(elem_file);
        if(!os.good()) {
            std::cout << "no file at " << elem_file << std::endl;
            return false;
        }

        for(auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it) {
            const libMesh::Elem &e = **it;

            for(int i = 0; i < e.n_nodes(); ++i) {
                int32_t idx = e.node_id(i);
                os.write((const char *)&idx, sizeof(idx));
            }
        }

        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_coords(const Path &folder, const libMesh::MeshBase &mesh, const std::string &postfix)
    {
        std::string coord_file = folder / ("coord_" + postfix);

        std::ofstream os;

        auto nodes_end   = mesh.active_nodes_end();
        auto nodes_begin = mesh.active_nodes_begin();

        int spatial_dim = mesh.spatial_dimension();

        for(int d = 0; d < spatial_dim; ++d) {

            os.open(coord_file + std::to_string(d) + ".raw");
            if(!os.good()) { assert(false); return false; }

            for(auto it = nodes_begin; it != nodes_end; ++it) {
                auto &n = **it;

                float x = n(d);
                os.write((const char *)(&x), sizeof(x));
            }

            os.close();
        }

        return true;
    }

    bool DiegoMeshWriter::write(const Path &folder, const libMesh::MeshBase &mesh)
    {
        return write_headers(folder, mesh) && write_elements(folder, mesh) && write_coords(folder, mesh);
    }

}