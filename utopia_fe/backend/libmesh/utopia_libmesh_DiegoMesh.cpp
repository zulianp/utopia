#include "utopia_libmesh_DiegoMesh.hpp"

#include <fstream>
#include <numeric>

#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/mesh.h"
#include "libmesh/node.h"
#include "libmesh/utility.h"

namespace utopia {

    bool DiegoMeshWriter::write_coords(const Path &folder,
                                       const libMesh::MeshBase &mesh,
                                       const Selector &selector,
                                       const std::string &postfix) {
        return write_coords(folder, mesh, selector.node_index, postfix);
    }

    void DiegoMeshWriter::build_index(const libMesh::MeshBase &mesh,
                                      const std::set<int> &blocks,
                                      int &n_elements,
                                      std::vector<int> &node_index,
                                      int &n_nodes) {
        node_index.resize(mesh.n_local_nodes());
        std::fill(node_index.begin(), node_index.end(), 0);

        n_elements = 0;

        for (auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it) {
            auto &e = **it;

            if (blocks.count(e.subdomain_id()) > 0) {
                ++n_elements;

                for (int i = 0; i < e.n_nodes(); ++i) {
                    int32_t idx = e.node_id(i);
                    node_index[idx] = 1;
                }
            }
        }

        n_nodes = std::accumulate(node_index.begin(), node_index.end(), 0);

        auto nodes_end = mesh.active_nodes_end();
        auto nodes_begin = mesh.active_nodes_begin();

        int current_offset = 0;
        // for(int i = 0; i < n_nodes; ++i) {

        for (auto it = nodes_begin; it != nodes_end; ++it) {
            auto &n = **it;

            if (node_index[n.id()] > 0) {
                node_index[n.id()] = 1 + current_offset++;
            }
        }
    }

    bool DiegoMeshWriter::write_coords(const Path &folder,
                                       const libMesh::MeshBase &mesh,
                                       const std::vector<int> &node_index,
                                       const std::string &postfix) {
        std::string coord_file = folder / ("coord_" + postfix);

        std::ofstream os;

        auto nodes_end = mesh.active_nodes_end();
        auto nodes_begin = mesh.active_nodes_begin();

        int spatial_dim = mesh.spatial_dimension();

        for (int d = 0; d < spatial_dim; ++d) {
            os.open(coord_file + std::to_string(d) + ".raw");
            if (!os.good()) {
                assert(false);
                return false;
            }

            for (auto it = nodes_begin; it != nodes_end; ++it) {
                auto &n = **it;

                if (node_index[n.id()] > 0) {
                    float x = n(d);
                    os.write((const char *)(&x), sizeof(x));
                }
            }

            os.close();
        }

        return true;
    }

    bool DiegoMeshWriter::write(const Path &folder, const libMesh::MeshBase &mesh, const std::set<int> &blocks) {
        if (blocks.empty()) {
            return write(folder, mesh);
        }

        std::string node_count_path = folder / "node_count.txt";
        std::string element_count_path = folder / "element_count.txt";
        std::string format_file = folder / "format.txt";
        std::string elem_file = folder / "elem.raw";

        std::vector<int> node_index;
        int n_elements = 0, n_nodes = 0;

        build_index(mesh, blocks, n_elements, node_index, n_nodes);

        std::ofstream os;
        os.open(node_count_path.c_str());

        if (!os.good()) {
            std::cout << "no file at " << node_count_path << std::endl;
            return false;
        }

        os << n_nodes << "\n";
        os.close();

        os.open(element_count_path.c_str());
        if (!os.good()) {
            std::cout << "no file at " << element_count_path << std::endl;
            return false;
        }

        os << n_elements << "\n";
        os.close();

        write_coords(folder, mesh, node_index);

        os.open(elem_file.c_str());

        if (!os.good()) {
            return false;
        }

        for (auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it) {
            auto &e = **it;

            if (blocks.count(e.subdomain_id()) > 0) {
                for (int i = 0; i < e.n_nodes(); ++i) {
                    assert(node_index[e.node_id(i)] > 0);

                    int32_t idx = node_index[e.node_id(i)] - 1;
                    os.write((const char *)&idx, sizeof(idx));
                }
            }
        }

        os.close();

        return true;
    }

    void DiegoMeshWriter::read(Input &in) {}

    bool DiegoMeshWriter::write_field(const Path &folder,
                                      const std::string &field_name,
                                      const libMesh::MeshBase &mesh,
                                      const UVector &field,
                                      const int component,
                                      const int sys_num) {
        std::string path = folder / (field_name + ".raw");

        std::ofstream os;

        os.open(path);
        if (!os.good()) {
            std::cout << "no file at " << path << std::endl;
            return false;
        }

        Read<UVector> r_d(field);

        auto m_it = mesh.local_nodes_begin();
        auto m_end = mesh.local_nodes_end();

        for (; m_it != m_end; ++m_it) {
            const int dof_id = (*m_it)->dof_number(sys_num, component, 0);
            float v = field.get(dof_id);
            os.write((char *)&v, sizeof(v));
        }

        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_headers(const Path &folder, const libMesh::MeshBase &mesh) {
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

        if (!os.good()) {
            std::cout << "no file at " << node_count_path << std::endl;
            return false;
        }

        os << mesh.n_nodes() << "\n";
        os.close();

        os.open(element_count_path.c_str());
        if (!os.good()) {
            std::cout << "no file at " << element_count_path << std::endl;
            return false;
        }

        os << mesh.n_active_elem() << "\n";
        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_elements(const Path &folder, const libMesh::MeshBase &mesh) {
        std::string elem_file = folder / "elem.raw";
        std::ofstream os;
        os.open(elem_file);
        if (!os.good()) {
            std::cout << "no file at " << elem_file << std::endl;
            return false;
        }

        for (auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it) {
            const libMesh::Elem &e = **it;

            for (int i = 0; i < e.n_nodes(); ++i) {
                int32_t idx = e.node_id(i);
                os.write((const char *)&idx, sizeof(idx));
            }
        }

        os.close();
        return true;
    }

    bool DiegoMeshWriter::write_coords(const Path &folder, const libMesh::MeshBase &mesh, const std::string &postfix) {
        std::string coord_file = folder / ("coord_" + postfix);

        std::ofstream os;

        auto nodes_end = mesh.active_nodes_end();
        auto nodes_begin = mesh.active_nodes_begin();

        int spatial_dim = mesh.spatial_dimension();

        for (int d = 0; d < spatial_dim; ++d) {
            os.open(coord_file + std::to_string(d) + ".raw");
            if (!os.good()) {
                assert(false);
                return false;
            }

            for (auto it = nodes_begin; it != nodes_end; ++it) {
                auto &n = **it;

                float x = n(d);
                os.write((const char *)(&x), sizeof(x));
            }

            os.close();
        }

        return true;
    }

    bool DiegoMeshWriter::write(const Path &folder, const libMesh::MeshBase &mesh) {
        return write_headers(folder, mesh) && write_elements(folder, mesh) && write_coords(folder, mesh);
    }

}  // namespace utopia