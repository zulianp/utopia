#ifndef UTOPIA_NORMALIZE_ROWS_HPP
#define UTOPIA_NORMALIZE_ROWS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class NormalizeRows {
    public:
        using Vector   = utopia::Wrapper<UTOPIA_VECTOR(Matrix), 1>;
        using Scalar   = UTOPIA_SCALAR(Matrix);
        using SizeType = UTOPIA_SIZE_TYPE(Matrix);

        static void apply(Matrix &in_out)
        {
            Vector row_sum = sum(in_out, 1);

            each_transform(row_sum, row_sum, [](const SizeType i, const Scalar val) -> Scalar {
               UTOPIA_UNUSED(i);

                if(val != 0.) {
                    return 1./val;
                } else {
                    return val;
                }
            }); 


            Matrix s = diag(row_sum);
            in_out = s * in_out;
        }
    };

    template<class Matrix>
    void normalize_rows(Matrix &in_out)
    {
        NormalizeRows<Matrix>::apply(in_out);
    }

}

#endif //UTOPIA_NORMALIZE_ROWS_HPP
