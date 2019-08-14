#include "utopia_libmesh_FunctionSpace.hpp"
#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"
#include "utopia_LibMeshBackend.hpp"

namespace utopia {

    void write(const Path &path, LibMeshFunctionSpace &space, UVector &x)
    {
        utopia::convert(x, *space.equation_system().solution);
        space.equation_system().solution->close();

        if(space.mesh().comm().size() == 1) {
            libMesh::ExodusII_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
        } else {
            libMesh::Nemesis_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
        }
    }

    
    std::size_t max_nnz_x_row(const LibMeshFunctionSpace &space)
    {
        return max_nnz_x_row(space.dof_map());
    }

    std::size_t max_nnz_x_row(const libMesh::DofMap &dof_map)
    {

        std::size_t nnz = 0;

        if(!dof_map.get_n_nz().empty()) {
            nnz = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

        }

        if(!dof_map.get_n_oz().empty()) {
            nnz += *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
        }


        return nnz;
    }
}
