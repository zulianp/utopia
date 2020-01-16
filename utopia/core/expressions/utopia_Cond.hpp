#ifndef UTOPIA_COND_HPP
#define UTOPIA_COND_HPP

#include "utopia_Traits.hpp"

namespace utopia {

    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class Cond {
    public:
        using Scalar = UTOPIA_SCALAR(Matrix);

        static Scalar apply(const Matrix &) 
        {
             static_assert(Backend < HOMEMADE, "Cond not implemented for this backend");
        }
    };

    template<class Matrix>
    inline UTOPIA_SCALAR(Matrix) cond(const Matrix &mat)
    {
        return Cond<Matrix>::apply(mat);
    }
}

#endif //UTOPIA_COND_HPP
