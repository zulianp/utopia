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
    void set_zero_rows(Wrapper<Tensor, 2> &w,  const std::vector<typename utopia::Traits<Tensor>::SizeType> & index)
    {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().set_zero_rows(w.implementation(), index);
    }

    // not sure how to name this one 
    template<class MatTensor, class VectorTensor>
    void apply_BC_to_system(Wrapper<MatTensor, 2> &A, Wrapper<VectorTensor, 1> &x, Wrapper<VectorTensor, 1> &rhs,  const std::vector<typename utopia::Traits<VectorTensor>::SizeType> & index)
    {
        Backend<typename Traits<MatTensor>::Scalar, Traits<MatTensor>::Backend>::Instance().apply_BC_to_system(A.implementation(), x.implementation(), rhs.implementation(), index);
    }

}

#endif //utopia_TEMP_HPP
