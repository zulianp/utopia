#include "utopia_libmesh_FunctionSpace.hpp"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/vtk_io.h"

#include "utopia_LibMeshBackend.hpp"

namespace utopia {

    void write(const Path &path, LibMeshFunctionSpace &space, UVector &x) {
        std::cout << "Writing solution at " << path.to_string() << std::endl;
        space.equation_system().solution->add(-6.0);
        std::cout << "x_norm: " << x.norm2() << std::endl;

        utopia::convert(x, *space.equation_system().solution);
        space.equation_system().solution->close();

        std::cout << "lm_x_norm: " << space.equation_system().solution->l2_norm() << std::endl;

        if (path.extension() == "pvtu") {
            libMesh::VTKIO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());

        } else {
            // if (space.mesh().comm().size() == 1) {
            libMesh::ExodusII_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
            // } else {
            //     libMesh::Nemesis_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
            // }
        }
    }

    std::size_t max_nnz_x_row(const LibMeshFunctionSpace &space) { return max_nnz_x_row(space.dof_map()); }

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
