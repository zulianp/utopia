#ifndef UTOPIA_TRANSFER_UTILS_HPP
#define UTOPIA_TRANSFER_UTILS_HPP

#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_MaxRowNNZ.hpp"

namespace utopia {

    inline static void tensorize(const USparseMatrix &T_x, const SizeType n_var, USparseMatrix &T)
    {
        auto max_nnz = utopia::max_row_nnz(T_x);
        T = local_sparse(local_size(T_x), max_nnz);

        Write<USparseMatrix> w(T);
        each_read(T_x, [&](const SizeType i, const SizeType j, const double value) {
            for(SizeType k = 0; k < n_var; ++k) {
                T.set(i + k, j + k, value);
            }
        });
    }

    inline static void tensorize(const SizeType n_var, UVector &t)
    {
        ReadAndWrite<UVector> w(t);
        auto r = range(t);

        for(auto i = r.begin(); i < r.end(); i += n_var) {
            const auto value = t.get(i);

            for(SizeType k = 1; k < n_var; ++k) {
                t.set(i + k, value);
            }
        }
    }


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
