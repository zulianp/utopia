#ifndef UTOPIA_TRANSFER_UTILS_HPP
#define UTOPIA_TRANSFER_UTILS_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {

    bool assemble_interpolation(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const int n_var = 1);

    bool assemble_projection(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const bool use_biorth = false, const int n_var = 1);

    bool assemble_coupling(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B);

    //use different lagr mult space
    bool assemble_projection(
        LibMeshFunctionSpace &from,
        LibMeshFunctionSpace &to,
        LibMeshFunctionSpace &lagr,
        USparseMatrix &B, USparseMatrix &D);


    bool assemble_interpolation(
        const libMesh::MeshBase &mesh, 
        const libMesh::DofMap &dof_map_p1,
        const libMesh::DofMap &dof_map_p2,
        USparseMatrix &mat
    );

}

#endif //UTOPIA_TRANSFER_UTILS_HPP
