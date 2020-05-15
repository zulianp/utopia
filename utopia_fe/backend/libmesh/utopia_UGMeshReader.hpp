#ifndef UG_MESH_READER_HPP
#define UG_MESH_READER_HPP

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace libMesh {
    class MeshBase;
}

namespace utopia {

    class UGXMeshReader {
    public:
        typedef std::vector<std::string> StringVector;

        virtual bool read(const std::string &path, libMesh::MeshBase &mesh) {
            std::ifstream file(path.c_str(), std::ifstream::in);
            if (!file.good()) {
                std::cerr << "[Error] unable to open file: " << path << std::endl;
                file.close();
                return false;
            }

            if (!read(file, mesh)) {
                return false;
            }

            file.close();
            return true;
        }

        bool read(std::istream &is, libMesh::MeshBase &mesh);
    };
}  // namespace utopia

#endif  // UG_MESH_READER_HPP
