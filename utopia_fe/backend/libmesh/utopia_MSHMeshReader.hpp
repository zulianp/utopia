#ifndef MSH_MESH_READER_HPP
#define MSH_MESH_READER_HPP 

#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>

namespace libMesh {
	class MeshBase;
}

namespace utopia {

	class MeshReaderOpts {
	public:

		MeshReaderOpts()
		: append_mode(false), subdomain_id(1)
		{}

		bool append_mode;
		std::size_t subdomain_id;
	};

	class MSHMeshReader {
	public:
		typedef std::vector<std::string> StringVector;

		
		bool read(std::istream &is, libMesh::MeshBase &mesh, const MeshReaderOpts &opts = MeshReaderOpts());

		bool read(const std::string &path, libMesh::MeshBase &mesh, const MeshReaderOpts &opts = MeshReaderOpts())
		{
			std::ifstream file(path.c_str(), std::ifstream::in);
			if(!file.good()) {
				std::cerr << "[Error] unable to open file: " << path << std::endl;
				file.close();
				return false;
			}

			if(!read(file, mesh, opts)) {
				file.close();
				return false;
			}

			file.close();
			return true;
		}

	};
}

#endif //MSH_MESH_READER_HPP

