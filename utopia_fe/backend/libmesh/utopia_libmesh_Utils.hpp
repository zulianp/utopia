#ifndef UTOPIA_LIBMESH_UTILS_HPP
#define UTOPIA_LIBMESH_UTILS_HPP

#include "libmesh/enum_elem_type.h"
#include "libmesh/libmesh_version.h"

#include <cassert>

namespace libMesh {
    class DofMap;
}

namespace utopia {
    inline bool is_pyramid(const int type) {
        return type == static_cast<int>(libMesh::PYRAMID5) || type == static_cast<int>(libMesh::PYRAMID13) ||
               type == static_cast<int>(libMesh::PYRAMID14);
    }

    inline bool is_prism(const int type) {
        return type == static_cast<int>(libMesh::PRISM6) || type == static_cast<int>(libMesh::PRISM15) ||
               type == static_cast<int>(libMesh::PRISM18);
    }

    inline bool is_hex(const int type) {
        return type == static_cast<int>(libMesh::HEX8) || type == static_cast<int>(libMesh::HEX20) ||
               type == static_cast<int>(libMesh::HEX27);
    }

    inline bool is_quad(const int type) {
        return type == static_cast<int>(libMesh::QUAD4) || type == static_cast<int>(libMesh::QUAD8) ||
               type == static_cast<int>(libMesh::QUAD9) || type == static_cast<int>(libMesh::QUADSHELL4) ||
               type == static_cast<int>(libMesh::QUADSHELL8);
    }

    inline bool is_tri(const int type) {
        return type == static_cast<int>(libMesh::TRI3) || type == static_cast<int>(libMesh::TRI6) ||
               type == static_cast<int>(libMesh::TRISHELL3);
    }

    inline bool is_tet(const int type) {
        return type == static_cast<int>(libMesh::TET4) || type == static_cast<int>(libMesh::TET10);
    }

    inline bool is_edge(const int type) {
        return type == static_cast<int>(libMesh::EDGE2) || type == static_cast<int>(libMesh::EDGE3) ||
               type == static_cast<int>(libMesh::EDGE4);
    }

    inline bool is_valid_elem_type(const int type) { return type < static_cast<int>(libMesh::INVALID_ELEM); }

    inline libMesh::ElemType side_type(const libMesh::ElemType &type) {
        using namespace libMesh;

        switch (type) {
            case QUAD4:
                return EDGE2;
            case QUAD8:
                return EDGE3;
            case QUAD9:
                return EDGE3;
            case QUADSHELL4:
                return EDGE2;
            case QUADSHELL8:
                return EDGE3;
            case TRI3:
                return EDGE2;
            case TRI6:
                return EDGE3;
            case TET4:
                return TRI3;
            case TET10:
                return TRI6;
            case HEX8:
                return QUAD4;
            case HEX20:
                return QUAD8;
            case HEX27:
                return QUAD9;
            case EDGE2:
                return NODEELEM;
            case EDGE3:
                return NODEELEM;
            case EDGE4:
                return NODEELEM;
            case PRISM6:
                return QUAD4;  // FIXME this is not always the case
            case PYRAMID5:
                return QUAD4;  // FIXME this is not always the case
            default: {
                assert(false && "add special case");
                return libMesh::INVALID_ELEM;
            }
        }
    }

    std::size_t max_nnz_x_row(const libMesh::DofMap &dof_map);
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_UTILS_HPP
