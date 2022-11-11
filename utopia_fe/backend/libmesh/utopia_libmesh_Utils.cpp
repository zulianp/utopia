#include "utopia_libmesh_Utils.hpp"

#include <libmesh/dof_map.h>

namespace utopia {
    std::size_t max_nnz_x_row(const libMesh::DofMap &dof_map) {
        std::size_t nnz = 0;

        if (!dof_map.get_n_nz().empty()) {
            nnz = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());
        }

        if (!dof_map.get_n_oz().empty()) {
            nnz += *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
        }

        return nnz;
    }
}  // namespace utopia
