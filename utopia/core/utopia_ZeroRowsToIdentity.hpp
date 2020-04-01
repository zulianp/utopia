#ifndef UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP
#define UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP

#include "utopia_Traits.hpp"
#include "utopia_Size.hpp"
#include "utopia_Temp.hpp"
#include "utopia_Algorithms.hpp"

#include <vector>
#include <cmath>


namespace utopia {

    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class ZeroRowsToIdentity {
    public:
        using Traits   = utopia::Traits<Matrix>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        static void apply(Matrix &A, const Scalar tol = 0.)
        {
            Size s = local_size(A);
            const std::size_t n = s.get(0);
            std::vector<bool> is_zero_row(n, true);
            const std::size_t r_begin = row_range(A).begin();

            each_read(A,[&is_zero_row, r_begin, tol](const SizeType i, const SizeType, const double val) {
                if(device::abs(val) > tol) {
                    is_zero_row[i - r_begin] = false;
                }
            });

            std::vector<SizeType> idx;
            idx.reserve(is_zero_row.size());

            for(std::size_t i = 0; i < n; ++i) {
                if(is_zero_row[i]) {
                    idx.push_back(i + r_begin);
                }
            }

            //hack in case some rows do not exist (relevant for petsc)
            if(Traits::Backend == PETSC && A.is_sparse()) {
                Matrix temp;
                temp.identity(layout(A), 0.0);
                A += temp;
            }

            set_zero_rows(A, idx, 1.);
        }
    };

    template<class Matrix>
    void zero_rows_to_identity(Matrix &A, const UTOPIA_SCALAR(Matrix) tol = 0.)
    {
        ZeroRowsToIdentity<Matrix>::apply(A, tol);
    }
}

#endif //UTOPIA_ZERO_ROWS_TO_IDENTITY_HPP
