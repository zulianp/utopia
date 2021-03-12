#include "utopia_MSHMeshReader.hpp"
// #include "rapidxml.hpp"
// #include "rapidxml_print.hpp"
// #include "rapidxml_utils.hpp"

#include <libmesh/const_function.h>
#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/petsc_vector.h>

#include <algorithm>

namespace utopia {

    static bool read_point(const std::string &string, std::size_t &index, libMesh::Point &xyz) {
        typedef libMesh::Real Real;
        return 4 == std::sscanf(string.c_str(), "%ld %lg %lg %lg", &index, &xyz(0), &xyz(1), &xyz(2));
    }

    static bool read_triangle(const std::string &node_string, libMesh::MeshBase &mesh) { return false; }

    static bool read_element(const std::string &node_string,
                             const std::size_t node_id_offset,
                             const MeshReaderOpts &opts,
                             libMesh::MeshBase &mesh) {
        std::istringstream ss(node_string);

        std::size_t elem_number = 0, elem_type = 0, n_tags = 0, n_nodes = 0;
        ss >> elem_number;
        ss >> elem_type;
        ss >> n_tags;

        std::unique_ptr<libMesh::Elem> elem;

        switch (elem_type) {
            case 4: {
                n_nodes = 4;
                elem = libMesh::Elem::build(libMesh::TET4);
                break;
            }
            default: {
                std::cerr << "[Error] Not implemented for elem_type = " << elem_type << std::endl;
                return false;
            }
        }

        if (!ss.good()) {
            std::cerr << "[Error] unable to read basic meta description for element " << elem_number << std::endl;
            return false;
        }

        std::vector<std::size_t> tags, nodes;
        tags.reserve(n_tags);
        nodes.reserve(n_nodes);

        std::size_t val = 0;
        for (std::size_t i = 0; i < n_tags && ss.good(); ++i) {
            ss >> val;
            tags.push_back(val);
        }

        for (std::size_t i = 0; i < n_nodes && ss.good(); ++i) {
            ss >> val;
            nodes.push_back(val);
        }

        if (n_nodes != nodes.size()) {
            std::cerr << "inconsistent format for element " << elem_number << std::endl;
        }

        for (std::size_t i = 0; i < n_nodes; ++i) {
            elem->set_node(i) = mesh.node_ptr(node_id_offset + nodes[i]);
        }

        elem->subdomain_id() = opts.subdomain_id;
        mesh.add_elem(elem.release());
        return true;
    }

    static const bool is_empty_line(const std::string &line) {
        return std::find_if_not(begin(line), end(line), [](unsigned char c) { return std::isspace(c); }) == line.end();
    }

    bool MSHMeshReader::read(std::istream &is, libMesh::MeshBase &mesh, const MeshReaderOpts &opts) {
        static const std::string nodes_begin_marker = "$Nodes";
        static const auto nodes_begin_marker_size = nodes_begin_marker.size();

        static const std::string nodes_end_marker = "$EndNodes";
        static const auto nodes_end_marker_size = nodes_end_marker.size();

        static const std::string elements_begin_marker = "$Elements";
        static const auto elements_begin_marker_size = elements_begin_marker.size();

        static const std::string elements_end_marker = "$EndElements";
        static const auto elements_end_marker_size = elements_end_marker.size();

        std::size_t node_id_offset = 0;
        if (!opts.append_mode) {
            mesh.clear();
        } else {
            node_id_offset = mesh.n_nodes();
        }

        std::string line;
        std::size_t n_nodes = 0;
        std::size_t n_elements = 0;

        while (is.good()) {
            std::getline(is, line);

            // allow empty newlines
            if (line.empty()) continue;
            if (is_empty_line(line)) continue;

            // begin read nodes
            if (line.compare(0, nodes_begin_marker_size, nodes_begin_marker) == 0) {
                if (!is.good()) {
                    std::cerr << "[Error] Bad file content. Maybe corrupt." << std::endl;
                    return false;
                }

                // std::getline(is, line);

                is >> n_nodes;

                std::size_t node_id = 0;

                libMesh::Point xyz;

                while (is.good()) {
                    std::getline(is, line);

                    // allow empty newlines
                    if (line.empty()) continue;
                    if (is_empty_line(line)) continue;

                    if (line.compare(0, nodes_end_marker_size, nodes_end_marker) == 0) {
                        // finished reading nodes
                        break;
                    }

                    if (!read_point(line, node_id, xyz)) {
                        std::cerr << "[Error] Bad node format. Line: " << std::endl;
                        std::cerr << line << std::endl;
                        return false;
                        continue;
                    }

                    mesh.add_point(xyz, node_id_offset + node_id);
                }
            }
            // end read nodes

            // begin read elements
            if (line.compare(0, elements_begin_marker_size, elements_begin_marker) == 0) {
                if (!is.good()) {
                    std::cerr << "[Error] Bad file content. Maybe corrupt." << std::endl;
                    return false;
                }

                is >> n_elements;
                std::size_t element_id = 0;

                while (is.good()) {
                    std::getline(is, line);

                    // allow empty newlines
                    if (line.empty()) continue;
                    if (is_empty_line(line)) continue;

                    if (line.compare(0, elements_end_marker_size, elements_end_marker) == 0) {
                        // finished reading elements
                        break;
                    }

                    if (!read_element(line, node_id_offset, opts, mesh)) {
                        std::cerr << "[Error] Bad element format. Line: " << std::endl;
                        std::cerr << line << std::endl;
                        return false;
                    }
                }
            }
            // end read elements
        }

        if (!opts.append_mode) {
            mesh.prepare_for_use();
        }

        return true;
    }
}  // namespace utopia
