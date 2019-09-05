/*
* @Author: Alena Kopanicakova
* @Date:   2016-09-12
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2016-09-14
*/

#ifndef utopia_TEMP_HPP
#define utopia_TEMP_HPP


#include "utopia_Expression.hpp"
#include "utopia_Evaluator.hpp"
#include "utopia_Assign.hpp"
#include "utopia_Traits.hpp"
#include "utopia_InPlace.hpp"
#include "utopia_Mutable.hpp"
#include "utopia_Readable.hpp"
#include "utopia_Writable.hpp"
#include "utopia_Ranged.hpp"

#include <iostream>
#include <type_traits>


namespace utopia
{


    template<class Tensor>
    void set_zero_rows(Wrapper<Tensor, 2> &w,  const std::vector<typename utopia::Traits<Tensor>::SizeType> & index, const double diag)
    {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().set_zero_rows(w.implementation(), index, diag);
    }

    // not sure how to name this one
    template<class MatTensor, class VectorTensor>
    void apply_BC_to_system(Wrapper<MatTensor, 2> &A, Wrapper<VectorTensor, 1> &x, Wrapper<VectorTensor, 1> &rhs,  const std::vector<typename utopia::Traits<VectorTensor>::SizeType> & index)
    {
        Backend<typename Traits<MatTensor>::Scalar, Traits<MatTensor>::Backend>::Instance().apply_BC_to_system(A.implementation(), x.implementation(), rhs.implementation(), index);
    }

    template<class Matrix, class Vector>
    void set_zero_rows(Wrapper<Matrix, 2> &w, const Wrapper<Vector, 1> &indicator, const double diag = 0.)
    {
        using VectorT  = utopia::Wrapper<Vector, 1>;
        using Scalar   = UTOPIA_SCALAR(VectorT);
        using SizeType = UTOPIA_SIZE_TYPE(VectorT);

        std::vector<SizeType> index;
        //index.reserve(local_size(indicator).get(0));

        each_read(indicator, [&index](const SizeType i, const Scalar value) {
            if(value == 1.) {
                index.push_back(i);
            }
        });

        set_zero_rows(w, index, diag);
    }


    template<class Vector, int Backend = Traits<Vector>::Backend>
    class EvalVecUniqueSortSerial
    {
        public:
            static void apply(const Wrapper<Vector, 1> &/*x*/, Wrapper<Vector, 1> & /*sorted */, const int /*used_values */ )
            {
                static_assert(Traits<Vector>::Backend==PETSC, "EvalVecUniqueSortSerial implemented just for petsc backend.");
            }
    };

    template<class Vector>
    void vec_unique_sort_serial(const Wrapper<Vector, 1> &x, Wrapper<Vector, 1> &sorted, const int used_values = -1)
    {
        EvalVecUniqueSortSerial<Vector>::apply(x, sorted, used_values);
    }


    template<class Matrix>
    void chop_abs(Wrapper<Matrix, 2> &A, const double eps)
    {
        Backend<typename Traits<Matrix>::Scalar, Traits<Matrix>::Backend>::Instance().chop(A.implementation(), eps);
    }


    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ChopSmallerThan
    {
        public:
            static void apply(const Wrapper<Matrix, 2> &/*x*/, const double & /*eps*/)
            {
                static_assert(Traits<Matrix>::Backend==PETSC, "ChopSmallerThan implemented just for petsc backend.");
            }
    };


    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ChopBiggerThan
    {
        public:
            static void apply(const Wrapper<Matrix, 2> & /*x*/, const double & /*eps*/)
            {
                static_assert(Traits<Matrix>::Backend==PETSC, "ChopBiggerThan implemented just for petsc backend.");
            }
    };


    template<class Matrix>
    void chop_smaller_than(Wrapper<Matrix, 2> &A, const double eps)
    {
        ChopSmallerThan<Matrix>::apply(A, eps);
    }    


    template<class Matrix>
    void chop_bigger_than(Wrapper<Matrix, 2> &A, const double eps)
    {
        ChopBiggerThan<Matrix>::apply(A, eps);
    }    


}

#endif //utopia_TEMP_HPP
