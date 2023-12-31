#ifndef UTOPIA_LIBMESH_TYPES_HPP
#define UTOPIA_LIBMESH_TYPES_HPP

#include <iostream>

#include "utopia_fe_kokkos_fix.hpp"

#include "utopia.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"

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

    // template<typename T>
    // class LibMeshAlgebraTraits {
    // public:
    //     typedef T Scalar;
    //     typedef libMesh::dof_id_type SizeType;

    //     using Matrix = utopia::LMDenseMatrix;
    //     using Vector = utopia::LMDenseVector;

    //     enum {
    //         Backend = LIBMESH_TAG
    //     };

    // };

    // UTOPIA_MAKE_TRAITS_DENSE_TPL_1(libMesh::DenseMatrix, LibMeshAlgebraTraits);
    // UTOPIA_MAKE_TRAITS_DENSE_TPL_1(libMesh::DenseVector, LibMeshAlgebraTraits);

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

        // IndexSet u_row_dofs(row_dofs.size()), u_col_dofs(col_dofs.size());

        // std::copy(row_dofs.begin(), row_dofs.end(), u_row_dofs.begin());
        // std::copy(col_dofs.begin(), col_dofs.end(), u_col_dofs.begin());

        // mat.add_matrix(u_row_dofs, u_col_dofs, block.entries());
    }

    template <typename Dofs>
    inline static void set_matrix(const libMesh::DenseMatrix<double> &block,
                                  const Dofs &row_dofs,
                                  const Dofs &col_dofs,
                                  USparseMatrix &mat) {
        // using IndexSet = Traits<UVector>::IndexSet;
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

        // IndexSet u_row_dofs(row_dofs.size()), u_col_dofs(col_dofs.size());

        // std::copy(row_dofs.begin(), row_dofs.end(), u_row_dofs.begin());
        // std::copy(col_dofs.begin(), col_dofs.end(), u_col_dofs.begin());

        // mat.set_matrix(u_row_dofs, u_col_dofs, block.entries());
    }

    template <typename Dofs>
    inline static void set_matrix(const LMDenseMatrix &block,
                                  const Dofs &row_dofs,
                                  const Dofs &col_dofs,
                                  USparseMatrix &mat) {
        // using IndexSet = Traits<UVector>::IndexSet;
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

        // IndexSet u_row_dofs(row_dofs.size()), u_col_dofs(col_dofs.size());

        // std::copy(row_dofs.begin(), row_dofs.end(), u_row_dofs.begin());
        // std::copy(col_dofs.begin(), col_dofs.end(), u_col_dofs.begin());

        // mat.set_matrix(u_row_dofs, u_col_dofs, block.entries());
    }

    inline static void add_vector(const LMDenseVector &block,
                                  const std::vector<libMesh::dof_id_type> &dofs,
                                  UVector &vec) {
        assert(block.size() == dofs.size());

        using IndexSet = Traits<UVector>::IndexSet;
        // for(uint i = 0; i < dofs.size(); ++i) {
        // 	vec.add(dofs[i], block(i));
        // }
        IndexSet index;
        index.insert(index.begin(), dofs.begin(), dofs.end());
        vec.add_vector(index, block.entries());
    }

    // template<typename Scalar>
    // class Backend<Scalar, LIBMESH_TAG> //: public ScalarBackend<Scalar>
    // {
    // public:
    //     typedef libMesh::dof_id_type SizeType;

    //     inline static Backend &Instance()
    //     {
    //         static Backend instance_;
    //         return instance_;
    //     }

    //     template<typename T>
    //     inline static void write_lock(T &, WriteMode) {}

    //     template<typename T>
    //     inline static void write_unlock(T &, WriteMode) {}

    //     template<typename T>
    //     inline static void read_lock(T &) {}

    //     template<typename T>
    //     inline static void read_unlock(T &) {}

    //     template<class Tensor>
    //     inline static void set(Tensor &t, const SizeType row, const SizeType col, const Scalar value)
    //     {
    //         t(row, col) = value;
    //     }

    //     template<class Tensor>
    //     inline static void add(Tensor &t, const SizeType row, const SizeType col, const Scalar value)
    //     {
    //         t(row, col) += value;
    //     }

    //     template<class Tensor>
    //     inline static Scalar get(Tensor &t, const SizeType row, const SizeType col)
    //     {
    //         return t(row, col);
    //     }

    //     template<class Tensor>
    //     inline static void set(Tensor &t, const SizeType row, const Scalar value)
    //     {
    //         t(row) = value;
    //     }

    //     template<class Tensor>
    //     inline static void add(Tensor &t, const SizeType row, const Scalar value)
    //     {
    //         t(row) += value;
    //     }

    //     template<class Tensor>
    //     inline static Scalar get(Tensor &t, const SizeType row)
    //     {
    //         return t(row);
    //     }

    //     template<typename T>
    //     inline static void assign_transposed(libMesh::DenseMatrix<T> &left, const libMesh::DenseMatrix<T> &right)
    //     {
    //         right.get_transpose(left);
    //     }

    //     template<typename T>
    //     inline static void build(libMesh::DenseMatrix<T> &m, const Size &size, const Zeros & /*values*/) {
    //         m.resize(size.get(0), size.get(1));
    //         m.zero();
    //     }

    //     template<typename T>
    //     inline static Range range(const libMesh::DenseVector<T> &v) {
    //         return Range(0, v.size());
    //     }

    //     template<typename T>
    //     inline static Range row_range(const libMesh::DenseMatrix<T> &m) {
    //         return Range(0, m.m());
    //     }

    //     template<typename T>
    //     inline static Range col_range(const libMesh::DenseMatrix<T> &m) {
    //         return Range(0, m.n());
    //     }

    //     template<typename T>
    //     inline static void resize(libMesh::DenseMatrix<T> &mat, const Size &size)
    //     {
    //         if(mat.m() != size.get(0) && mat.n() != size.get(1)) {
    //             mat.resize(size.get(0), size.get(1));
    //         }
    //     }

    //     template<typename T>
    //     inline static void resize(libMesh::DenseVector<T> &vec, const Size &size)
    //     {
    //         if(vec.size() != size.get(0)) {
    //             vec.resize(size.get(0));
    //         }
    //     }

    //     template<class L, class R>
    //     inline static void assign(L &left, R &&right) {
    //         left = std::forward<R>(right);
    //     }

    //     template<class Tensor>
    //     inline static double norm2(const Tensor &)
    //     {
    //         assert(false);
    //         return 1.;
    //     }

    //     inline static double norm2(const LMDenseVector &vec)
    //     {
    //         return vec.l2_norm();
    //     }

    // private:
    //     Backend() {}

    // };

    // template<typename Result>
    // class Evaluator<Result, LIBMESH_TAG> {
    // public:

    //     template<class Left, class Right>
    //     static void eval(const Construct<Left, Right> &expr)
    //     {
    //         ExprInliner<Construct<Left, Right>>::eval(expr);
    //     }

    //     template<class Left, class Right>
    //     static void eval(const Assign<Left, Right> &expr)
    //     {
    //         ExprInliner<Assign<Left, Right>>::eval(expr);
    //     }

    //     template<class Left, class Right, class Operation>
    //     static void eval(const InPlace<Left, Right, Operation> &expr)
    //     {
    //         ExprInliner<InPlace<Left, Right, Operation>>::eval(expr);
    //     }

    //     template<class Derived>
    //     static auto eval(const Expression<Derived> &expr) -> decltype( ExprInliner<Derived>::eval(expr.derived()) )
    //     {
    //         return ExprInliner<Derived>::eval(expr.derived());
    //     }

    //     template<class Derived>
    //     inline static void eval(const Tensor<Derived, 2> &t, Size &size)
    //     {
    //         size.set_dims(2);
    //         size.set(0, t.derived().m());
    //         size.set(1, t.derived().n());
    //     }

    //     template<class Derived>
    //     inline static void eval(const Tensor<Derived, 1> &t, Size &size)
    //     {
    //         size.set_dims(1);
    //         size.set(0, t.derived().size());
    //     }

    //     template<class Derived>
    //     inline static libMesh::Real eval(const Norm<Tensor<Derived, 1>,2> &t)
    //     {
    //         return Backend<libMesh::Real, LIBMESH_TAG>::norm2(t.expr().derived());
    //     }
    // };

    // inline void disp(const LMDenseMatrix &mat, std::ostream &os = std::cout)
    // {
    //     mat.implementation().print(os);
    // }

    // inline void disp(const LMDenseVector &vec, std::ostream &os = std::cout)
    // {
    //     vec.implementation().print(os);
    // }

    // inline const libMesh::DenseMatrix<libMesh::Real> &raw_type(const Wrapper<libMesh::DenseMatrix<libMesh::Real>, 2>
    // &utopiaType)
    // {
    //     return utopiaType.implementation();
    // }

    // inline libMesh::DenseMatrix<libMesh::Real> &raw_type(Wrapper<libMesh::DenseMatrix<libMesh::Real>, 2> &utopiaType)
    // {
    //     return utopiaType.implementation();
    // }
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_TYPES_HPP
