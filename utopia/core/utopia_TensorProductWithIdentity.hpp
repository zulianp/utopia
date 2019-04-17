#ifndef UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP
#define UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Factory.hpp"

namespace utopia {
    template<class Matrix, int Backend = Traits<Matrix>::Backend>
    class TensorProductWithIdentity {
    public:
        using SizeType = UTOPIA_SIZE_TYPE(Matrix);
        using Scalar   = UTOPIA_SCALAR(Matrix);

        static void apply(const Matrix &in, const SizeType dim, Matrix &out)
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

            out = local_sparse(ls.get(0) * dim, ls.get(1) * dim, max_nnz_x_row);

            Write<Matrix> w(out);

            each_read(in, [&](const SizeType i, const SizeType j, const Scalar val) {
                for(SizeType d = 0; d < dim; ++d) {
                    out.set(i*dim + d, j*dim + d, val);
                }
            });
        }
    };

    template<class Matrix>
    void tensor_prod_with_identity(const Matrix &in, const UTOPIA_SIZE_TYPE(Matrix) dim, Matrix &out)
    {
        TensorProductWithIdentity<Matrix>::apply(in, dim, out);
    }
}

#endif //UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP
