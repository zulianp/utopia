#ifndef UTOPIA_MAX_ROW_NNZ_HPP
#define UTOPIA_MAX_ROW_NNZ_HPP

//#include "utopia_Each.hpp"
#include "utopia_Traits.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

namespace utopia {

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class MaxRowNNZ {
    public:
        using SizeType = typename Traits<Matrix>::SizeType;
        using Scalar = typename Traits<Matrix>::Scalar;

        static SizeType apply(const Matrix &in) {
            auto rr = row_range(in);
            auto ls = local_size(in);
            auto gs = size(in);

            std::vector<SizeType> nnz(rr.extent(), 0);

            in.derived().read([&](const SizeType i, const SizeType &, const Scalar &) { ++nnz[i - rr.begin()]; });

            SizeType max_nnz_x_row = 0;

            if (!nnz.empty()) {
                max_nnz_x_row = *std::max_element(nnz.begin(), nnz.end());
            }

            return max_nnz_x_row;
        }
    };

    template <class Derived>
    inline typename Traits<Derived>::SizeType max_row_nnz(const Tensor<Derived, 2> &mat) {
        return MaxRowNNZ<Tensor<Derived, 2>>::apply(mat.derived());
    }

}  // namespace utopia

#endif  // UTOPIA_MAX_ROW_NNZ_HPP
