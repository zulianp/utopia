#ifndef UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP
#define UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP

#include "utopia_Traits.hpp"
#include "utopia_Size.hpp"

#include <vector>

namespace utopia {

    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ZeroRowsToIdentity {
    public:
        static void apply(Matrix &A)
        {
            Size s = local_size(A);
            const std::size_t n = s.get(0);
            std::vector<bool> d(n, false);
            const std::size_t r_begin = row_range(A).begin();

            each_read(A,[&d, r_begin](const SizeType i, const SizeType, const double val) {
                if(val != 0.) {
                    d[i - r_begin] = true;
                }
            });

            Write<Matrix> w(A);
            for(std::size_t i = 0; i < n; ++i) {
                auto i_global = i + r_begin;
                A.set(i_global, i_global, 1.);
            }
        }
    };

    template<class Matrix>
    void zero_rows_to_identity(Matrix &A)
    {
        ZeroRowsToIdentity<Matrix>::apply(A);
    }
}

#endif //UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP
