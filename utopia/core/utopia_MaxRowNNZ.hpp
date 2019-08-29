#ifndef UTOPIA_MAX_ROW_NNZ_HPP
#define UTOPIA_MAX_ROW_NNZ_HPP

#include "utopia_Traits.hpp"
#include "utopia_Each.hpp"
#include <vector>

namespace utopia {

    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class MaxRowNNZ {
    public:
        using SizeType = UTOPIA_SIZE_TYPE(Matrix);
        using Scalar   = UTOPIA_SCALAR(Matrix);

        static SizeType apply(const Matrix &in)
        {
            auto rr = row_range(in);
            auto ls = local_size(in);
            auto gs = size(in);

            std::vector<SizeType> nnz(rr.extent(), 0);
           
            each_read(in, [&](const SizeType i, const SizeType j, const Scalar val) {
                UTOPIA_UNUSED(j);
                UTOPIA_UNUSED(val);

                ++nnz[i - rr.begin()];
            });

            SizeType max_nnz_x_row = 0;

            if(!nnz.empty()) {
                max_nnz_x_row = *std::max_element(nnz.begin(), nnz.end());
            }

            return max_nnz_x_row;
        }
    };

    template<class Tensor> 
    inline SizeType max_row_nnz(const Wrapper<Tensor, 2> &mat)
    {
        return MaxRowNNZ<Wrapper<Tensor, 2>>::apply(mat);
    }

}

#endif //UTOPIA_MAX_ROW_NNZ_HPP
