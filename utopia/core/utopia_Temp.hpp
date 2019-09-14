#ifndef UTOPIA_TEMP_HPP
#define UTOPIA_TEMP_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Expression.hpp"
// #include "utopia_Evaluator.hpp"
// #include "utopia_Assign.hpp"
// #include "utopia_Traits.hpp"
// #include "utopia_InPlace.hpp"
// #include "utopia_Mutable.hpp"
// #include "utopia_Readable.hpp"
// #include "utopia_Writable.hpp"
// #include "utopia_Ranged.hpp"

#include <iostream>
#include <type_traits>

namespace utopia {

    template<class Derived>
    void set_zero_rows(
        Tensor<Derived, 2> &w, 
        const typename Traits<Derived>::IndexSet &index,
        const typename Traits<Derived>::Scalar diag)
    {
        w.derived().set_zero_rows(index, diag);
    }


    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class ApplyBCToSystem {
    public:
        template<class IndexSetT>
        void apply(Matrix &, Vector &, Vector&, const IndexSetT &)
        {
            static_assert(Backend < HOMEMADE, "implement for specific backend");
        }
    };

    //FIXME not sure how to name this one: change name to apply_equality_constraints_to_system
    template<class MatDerived, class VecDerived>
    void apply_BC_to_system(
        Tensor<MatDerived, 2> &A,
        Tensor<VecDerived, 1> &x,
        Tensor<VecDerived, 1> &rhs,
        const typename Traits<MatDerived>::IndexSet &constrained_idx)
    {
        ApplyBCToSystem<MatDerived, VecDerived>::apply(A.derived(), x.derived(), rhs.derived(), constrained_idx);
    }

    template<class Matrix, class Vector>
    void set_zero_rows(Tensor<Matrix, 2> &w, const Tensor<Vector, 1> &indicator, const double diag = 0.)
    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        //FIXME maybe use array type once is available
        using IndexSet = typename Traits<Vector>::IndexSet;

        IndexSet index;
        //index.reserve(local_size(indicator).get(0));

        {
            //Write<IndexSet> w(index);
            each_read(indicator.derived(), [&index](const SizeType i, const Scalar value) {
                if(value == 1.) {
                    index.push_back(i);
                }
            });
        }

        set_zero_rows(w, index, diag);
    }


    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalVecUniqueSortSerial
    {
        public:
            static void apply(const Tensor<Vector, 1> &/*x*/, Tensor<Vector, 1> & /*sorted */, const int /*used_values */ )
            {
                static_assert(Traits<Vector>::Backend==PETSC, "EvalVecUniqueSortSerial implemented just for petsc backend.");
            }
    };

    template<class Vector>
    void vec_unique_sort_serial(const Tensor<Vector, 1> &x, Tensor<Vector, 1> &sorted, const int used_values = -1)
    {
        EvalVecUniqueSortSerial<Vector>::apply(x, sorted, used_values);
    }


    template<class Matrix>
    void chop_abs(Tensor<Matrix, 2> &A, const double eps)
    {
        Backend<typename Traits<Matrix>::Scalar, Traits<Matrix>::Backend>::Instance().chop(A.implementation(), eps);
    }


    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ChopSmallerThan
    {
        public:
            static void apply(const Tensor<Matrix, 2> &/*x*/, const double & /*eps*/)
            {
                static_assert(Traits<Matrix>::Backend<HOMEMADE, "ChopSmallerThan implemented just for petsc backend.");
            }
    };


    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ChopGreaterThan
    {
        public:
            static void apply(const Tensor<Matrix, 2> & /*x*/, const double & /*eps*/)
            {
                static_assert(Traits<Matrix>::Backend<HOMEMADE, "ChopGreaterThan implemented just for petsc backend.");
            }
    };


    template<class Matrix>
    void chop_smaller_than(Tensor<Matrix, 2> &A, const double eps)
    {
        ChopSmallerThan<Matrix>::apply(A, eps);
    }    


    template<class Matrix>
    void chop_bigger_than(Tensor<Matrix, 2> &A, const double eps)
    {
        ChopGreaterThan<Matrix>::apply(A, eps);
    }    

}

#endif //UTOPIA_TEMP_HPP
