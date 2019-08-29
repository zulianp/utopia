#ifndef UTOPIA_MOONOLITH_PERMUTATIONS_HPP
#define UTOPIA_MOONOLITH_PERMUTATIONS_HPP

#include "utopia.hpp"
#include "utopia_fe_base.hpp"

#include "moonolith_mesh.hpp"
#include "moonolith_function_space.hpp"

#include <vector>

namespace utopia {
    template<int Dim>
    using MoonolithFunctionSpace = moonolith::FunctionSpace<moonolith::Mesh<double, Dim>>;
    
    template<int Dim>
    void make_permutation(
        const MoonolithFunctionSpace<Dim> &from, 
        const MoonolithFunctionSpace<Dim> &to, 
        USparseMatrix &mat)
    {
        using Scalar   = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar>   vals(1, 1.0);

        auto max_nnz = from.dof_map().max_nnz(); assert(max_nnz > 0);
        mat = local_sparse(to.dof_map().n_local_dofs(), from.dof_map().n_local_dofs(), max_nnz);

        std::size_t n_elems = from.dof_map().n_elements();

        assert(n_elems == to.dof_map().n_elements());

        Write<USparseMatrix>  w(mat, utopia::GLOBAL_INSERT);

        for(std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map().dofs(e);
            const auto &to_dofs   = to.dof_map().dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to   = to_dofs.size();

            assert(n_from == n_to);

            for(std::size_t i = 0; i < n_to; ++i) {
                irows[0] = to_dofs[i];
                icols[0] = from_dofs[i];
                mat.set_matrix(irows, icols, vals);
            }
        }
    }
    
    template<int Dim>
    void make_vector_permutation(
     const int dim,
     const MoonolithFunctionSpace<Dim> &from,
     const MoonolithFunctionSpace<Dim> &to,
     USparseMatrix &mat)
    {
        using Scalar   = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar>   vals(1, 1.0);

        auto max_nnz =  from.dof_map().max_nnz(); assert(max_nnz > 0);
        mat = local_sparse(to.dof_map().n_local_dofs(), from.dof_map().n_local_dofs() * dim, max_nnz);

        std::size_t n_elems = from.dof_map().n_elements();

        assert(n_elems == to.dof_map().n_elements());

        Write<USparseMatrix>  w(mat, utopia::GLOBAL_INSERT);

        for(std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map().dofs(e);
            const auto &to_dofs   = to.dof_map().dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to   = to_dofs.size();

            assert(n_from == n_to);

            for(std::size_t i = 0; i < n_to; ++i) {
                for(int d = 0; d < dim; ++d) {
                    irows[0] = to_dofs[i] + d;
                    icols[0] = from_dofs[i] * dim + d;
                    mat.set_matrix(irows, icols, vals);
                }
            }
        }
    }
    
    template<int Dim>
    void make_tensorize_permutation(
        const int dim,
        const MoonolithFunctionSpace<Dim> &from,
        const MoonolithFunctionSpace<Dim> &to,
        USparseMatrix &mat)
    {
        using Scalar   = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar>   vals(1, 1.0);

        auto max_nnz = from.dof_map().max_nnz(); assert(max_nnz > 0);
        mat = local_sparse(to.dof_map().n_local_dofs() * dim, from.dof_map().n_local_dofs(), max_nnz);

        std::size_t n_elems = from.dof_map().n_elements();

        assert(n_elems == to.dof_map().n_elements());

        Write<USparseMatrix>  w(mat, utopia::GLOBAL_INSERT);

        for(std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map().dofs(e);
            const auto &to_dofs   = to.dof_map().dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to   = to_dofs.size();

            assert(n_from == n_to);

            for(std::size_t i = 0; i < n_to; ++i) {
                for(int d = 0; d < dim; ++d) {
                    irows[0] = to_dofs[i] + d;
                    icols[0] = from_dofs[i];
                    mat.set_matrix(irows, icols, vals);
                }
            }
        }
    }
}

#endif //UTOPIA_MOONOLITH_PERMUTATIONS_HPP
