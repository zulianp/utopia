#ifndef UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP
#define UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP

#include "utopia_Base.hpp"
#include "utopia_Factory.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

#include <algorithm>

namespace utopia {
    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class TensorProductWithIdentity {
    public:
        using Traits = utopia::Traits<Matrix>;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using Comm = typename Traits::Communicator;

        static void apply(const Matrix &in, const SizeType dim, Matrix &out) {
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

            if (!nnz.empty()) {
                max_nnz_x_row = *std::max_element(nnz.begin(), nnz.end());
            }

            // FIXME
            out.sparse(layout(in.comm(), ls.get(0) * dim, ls.get(1) * dim, Traits::determine(), Traits::determine()),
                       max_nnz_x_row,
                       max_nnz_x_row);

            Write<Matrix> w(out);

            each_read(in, [&](const SizeType i, const SizeType j, const Scalar val) {
                const auto offset_i = i * dim;
                const auto offset_j = j * dim;

                for (SizeType d = 0; d < dim; ++d) {
                    out.set(offset_i + d, offset_j + d, val);
                }
            });
        }
    };

    template <class Matrix>
    void tensor_prod_with_identity(const Matrix &in, const UTOPIA_SIZE_TYPE(Matrix) dim, Matrix &out) {
        TensorProductWithIdentity<Matrix>::apply(in, dim, out);
    }
}  // namespace utopia

#endif  // UTOPIA_TENSOR_PRODUCT_WITH_IDENTITY_HPP
