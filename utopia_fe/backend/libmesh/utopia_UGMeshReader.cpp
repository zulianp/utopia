#include "utopia_UGMeshReader.hpp"
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"
// #include "rapidxml_iterators.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/elem.h>

namespace utopia {

    static bool read_nodes(const int n_dims, const std::string &node_string, libMesh::MeshBase &mesh)
    {
        typedef libMesh::Real Real;

        std::vector<Real> pnt;
        pnt.reserve(n_dims);

        std::size_t offset = 0;
        std::size_t index = 0;
        std::size_t vertex_id = 0;

        std::string num;
        bool stop = false;
        while(offset != std::string::npos && !stop) {
            pnt.clear();

            for(std::size_t i = 0; i < n_dims; ++i) {
                index = node_string.find_first_of(" \n\t", offset);
                num = node_string.substr(offset, index - offset);
                pnt.push_back(atof(num.c_str()));

                if(index == std::string::npos) {
                    stop = true;
                }

                offset = index + 1;
            }

            if(n_dims == 3) {
                libMesh::Point xyz(pnt[0], pnt[1], pnt[2]);
                mesh.add_point(xyz, vertex_id++);
                // std::cout << vertex_id << "\t"; xyz.print(std::cout); std::cout << std::endl;
            } else {
                assert(false);
            }
        }

        return true;
    }

    static bool read_triangles(const std::string &node_string, libMesh::MeshBase &mesh)
    {
        return false;
    }

    static bool read_tetrahedra(const std::string &node_string, libMesh::MeshBase &mesh)
    {
        typedef libMesh::Real Real;

        const std::size_t n_nodes = 4;

        std::size_t offset = 0;
        std::size_t index = 0;
        std::size_t element_id = 0;

        std::string num;
        bool stop = false;
        while(offset != std::string::npos && !stop) {
            auto elem = libMesh::Elem::build(libMesh::TET4).release();

            for(std::size_t i = 0; i < n_nodes; ++i) {
                index = node_string.find_first_of(" \n\t", offset);
                num = node_string.substr(offset, index - offset);

                elem->set_node(i) = mesh.node_ptr(atoi(num.c_str()));

                if(index == std::string::npos) {
                    stop = true;
                }

                offset = index + 1;
            }

            mesh.add_elem(elem);
        }

        return true;
    }

    bool UGXMeshReader::read(std::istream &is, libMesh::MeshBase &mesh)
    {
        using namespace rapidxml;

        file<> f(is);
        xml_document<> doc;
        doc.parse<0>(f.data());

        xml_node<> *grid_node = doc.first_node("grid");

        if(!grid_node) {
            std::cerr << "could not find <grid> node" << std::endl;
            return false;
        }

        xml_node<> *vertices_node = grid_node->first_node("vertices");
        if(!vertices_node) {
            std::cerr << "could not find <vertices> node" << std::endl;
            return false;
        }

        xml_attribute<> *dim_attr = vertices_node->first_attribute("coords");
        if(!dim_attr) {
            return false;
        }

        int n_dims = atoi(dim_attr->value());
        std::cout << "n_dims: " << n_dims << std::endl;

        if(!read_nodes(n_dims, vertices_node->value(), mesh)) {
            return false;
        }

        xml_node<> *tets_node = grid_node->first_node("tetrahedrons");
        if(!tets_node) {
            std::cerr << "could not find <tetrahedrons> node" << std::endl;
            return false;
        }

        if(!read_tetrahedra(tets_node->value(), mesh)) {
            return false;
        }

        mesh.prepare_for_use();
        return true;
    }
}

