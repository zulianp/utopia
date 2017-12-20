#ifndef UTOPIA_LIBMESH_UTILS_HPP
#define UTOPIA_LIBMESH_UTILS_HPP 

#include "libmesh/enum_elem_type.h"

namespace utopia {
	inline bool is_hex(const int type)
	{
		return type == static_cast<int>(libMesh::HEX8) || type == static_cast<int>(libMesh::HEX20) ||
		type == static_cast<int>(libMesh::HEX27);
	}
	
	inline bool is_quad(const int type)
	{
		return type == static_cast<int>(libMesh::QUAD4) || type == static_cast<int>(libMesh::QUADSHELL4) ||
		type == static_cast<int>(libMesh::QUAD8) || type == static_cast<int>(libMesh::QUAD9);
	}
	
	inline bool is_tri(const int type)
	{
		return type == static_cast<int>(libMesh::TRI3)  || type == static_cast<int>(libMesh::TRISHELL3) ||
		type == static_cast<int>(libMesh::TRI6);
	}
	
	inline bool is_tet(const int type)
	{
		return type == static_cast<int>(libMesh::TET4) || type == static_cast<int>(libMesh::TET10);
	}
	
	inline bool is_valid_elem_type(const int type)
	{
		return type < static_cast<int>(libMesh::INVALID_ELEM);
	}
}


#endif //UTOPIA_LIBMESH_UTILS_HPP
