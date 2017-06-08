#ifndef UG_MESH_READER_HPP
#define UG_MESH_READER_HPP 

#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>

namespace libMesh {
	class Mesh;
}

namespace utopia {

	class UGXMeshReader {
	public:
		typedef std::vector<std::string> StringVector;
		typedef libMesh::Mesh Mesh;
		

		virtual bool read(const std::string &path, Mesh &mesh) 
		{
			std::ifstream file(path.c_str(), std::ifstream::in);
			if(!file.good()) {
				std::cerr << "[Error] unable to open file: " << path << std::endl;
				file.close();
				return false;
			}

			if(!read(file, mesh)) {
				return false;
			}

			file.close();
			return true;
		}

		bool read(std::istream &is, Mesh &mesh);
	};
}

#endif //UG_MESH_READER_HPP

