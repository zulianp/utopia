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
            // does not link with new libmesh
            utopia_error("utopia_libmesh_FunctionSpace:: write function:: add support for VTK \n");
            // libMesh::VTKIO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());

        } else {
            // if (space.mesh().comm().size() == 1) {
            libMesh::ExodusII_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
            // } else {
            //     libMesh::Nemesis_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
            // }
        }
    }

    std::size_t max_nnz_x_row(const LibMeshFunctionSpace &space) { return max_nnz_x_row(space.dof_map()); }

}  // namespace utopia
