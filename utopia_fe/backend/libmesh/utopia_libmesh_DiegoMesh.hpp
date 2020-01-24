#ifndef UTOPIA_LIBMESH_DIEGOMESH_HPP
#define UTOPIA_LIBMESH_DIEGOMESH_HPP

#include "utopia_Path.hpp"
#include "libmesh/mesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/utility.h"

#include <fstream>

#include "utopia_fe_base.hpp"

namespace utopia {

    class DiegoMeshWriter {
    public:

        bool write_field(const Path &folder, const libMesh::MeshBase &mesh, const UVector &field)
        {
            return false;
        }

        bool write(const Path &folder, const libMesh::MeshBase &mesh)
        {
            //paths
            std::string node_count_path = folder / "node_count.txt";
            std::string element_count_path = folder / "element_count.txt";
            std::string coord_file = folder / "coord_";
            std::string elem_file  = folder / "elem.raw";
            std::string format_file = folder / "format.txt";

            std::ofstream os;

            os.open(format_file);
            if(!os.good()) {
                std::cout << "no file at " << format_file << std::endl;
                return false;
            }

            os << libMesh::Utility::enum_to_string((*mesh.active_local_elements_begin())->type()) << "\n";

            os.close();

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

    };
}

#endif //UTOPIA_LIBMESH_DIEGOMESH_HPP
