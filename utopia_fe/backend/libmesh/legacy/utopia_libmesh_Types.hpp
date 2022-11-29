#ifndef UTOPIA_LIBMESH_TYPES_HPP
#define UTOPIA_LIBMESH_TYPES_HPP

#include <iostream>

#include "utopia_fe_kokkos_fix.hpp"

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
// #include "utopia_libmesh_FEForwardDeclarations.hpp"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"

namespace utopia {

    typedef utopia::BlasMatrix<libMesh::Real> LMDenseMatrix;
    typedef utopia::BlasVector<libMesh::Real> LMDenseVector;

    template <class It>
    void print_vector(const It &begin, const It &end, std::ostream &os = std::cout) {
        for (It it = begin; it != end; ++it) {
            os << *it << " ";
        }

        os << "\n";
    }

    inline static void add_matrix(const libMesh::DenseMatrix<libMesh::Real> &block,
                                  const std::vector<libMesh::dof_id_type> &row_dofs,
                                  const std::vector<libMesh::dof_id_type> &col_dofs,
                                  libMesh::SparseMatrix<libMesh::Real> &mat) {
        mat.add_matrix(block, row_dofs, col_dofs);
    }

    inline static void add_vector(const libMesh::DenseVector<libMesh::Real> &block,
                                  const std::vector<libMesh::dof_id_type> &dofs,
                                  libMesh::NumericVector<libMesh::Real> &vec) {
        assert(block.size() == dofs.size());
        vec.add_vector(block, dofs);
    }

    inline static void get_vector(const UVector &vec,
                                  const std::vector<libMesh::dof_id_type> &dofs,
                                  libMesh::DenseVector<libMesh::Real> &el_vec) {
        el_vec.resize(dofs.size());
        int i = 0;
        for (auto test : dofs) {
            el_vec(i++) = vec.get(test);
        }
    }

    inline static void add_matrix(const LMDenseMatrix &block,
                                  const std::vector<libMesh::dof_id_type> &row_dofs,
                                  const std::vector<libMesh::dof_id_type> &col_dofs,
                                  USparseMatrix &mat) {
        using IndexSet = Traits<USparseMatrix>::IndexSet;
        Size s = size(mat);
        for (uint i = 0; i < row_dofs.size(); ++i) {
            for (uint j = 0; j < col_dofs.size(); ++j) {
                const libMesh::Real val = block.get(i, j);
                if (val != 0.0) {
                    assert(row_dofs[i] < s.get(0));
                    assert(col_dofs[j] < s.get(1));

                    mat.c_add(row_dofs[i], col_dofs[j], val);
                }
            }
        }
    }

    template <typename Dofs>
    inline static void set_matrix(const libMesh::DenseMatrix<double> &block,
                                  const Dofs &row_dofs,
                                  const Dofs &col_dofs,
                                  USparseMatrix &mat) {
        Size s = size(mat);
        for (uint i = 0; i < row_dofs.size(); ++i) {
            for (uint j = 0; j < col_dofs.size(); ++j) {
                const libMesh::Real val = block(i, j);
                if (val != 0.0) {
                    assert(row_dofs[i] < s.get(0));
                    assert(col_dofs[j] < s.get(1));

                    mat.c_set(row_dofs[i], col_dofs[j], val);
                }
            }
        }
    }

    template <typename Dofs>
    inline static void set_matrix(const LMDenseMatrix &block,
                                  const Dofs &row_dofs,
                                  const Dofs &col_dofs,
                                  USparseMatrix &mat) {
        Size s = size(mat);
        for (uint i = 0; i < row_dofs.size(); ++i) {
            for (uint j = 0; j < col_dofs.size(); ++j) {
                const libMesh::Real val = block.get(i, j);
                if (val != 0.0) {
                    assert(row_dofs[i] < s.get(0));
                    assert(col_dofs[j] < s.get(1));

                    mat.c_set(row_dofs[i], col_dofs[j], val);
                }
            }
        }
    }

    inline static void add_vector(const LMDenseVector &block,
                                  const std::vector<libMesh::dof_id_type> &dofs,
                                  UVector &vec) {
        assert(block.size() == dofs.size());

        using IndexSet = Traits<UVector>::IndexSet;

        IndexSet index;
        index.insert(index.begin(), dofs.begin(), dofs.end());
        vec.add_vector(index, block.entries());
    }

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_TYPES_HPP
