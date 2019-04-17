#ifndef UTOPIA_TRANSFER_UTILS_HPP
#define UTOPIA_TRANSFER_UTILS_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {

    bool assemble_interpolation(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D);

    bool assemble_projection(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B, USparseMatrix &D, const bool use_biorth = false);

    bool assemble_coupling(LibMeshFunctionSpace &from, LibMeshFunctionSpace &to, USparseMatrix &B);

    //use different lagr mult space
    bool assemble_projection(
        LibMeshFunctionSpace &from,
        LibMeshFunctionSpace &to,
        LibMeshFunctionSpace &lagr,
        USparseMatrix &B, USparseMatrix &D);

}

#endif //UTOPIA_TRANSFER_UTILS_HPP
