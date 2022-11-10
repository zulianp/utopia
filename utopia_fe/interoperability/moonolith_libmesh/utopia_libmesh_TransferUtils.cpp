#include "utopia_libmesh_TransferUtils.hpp"

#include "utopia_MaxRowNNZ.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    void tensorize(const USparseMatrix &T_x, const SizeType n_var, USparseMatrix &T) {
        auto max_nnz = utopia::max_row_nnz(T_x);
        T = local_sparse(local_size(T_x), max_nnz);

        Write<USparseMatrix> w(T);
        each_read(T_x, [&](const SizeType i, const SizeType j, const double value) {
            for (SizeType k = 0; k < n_var; ++k) {
                T.set(i + k, j + k, value);
            }
        });
    }

    void tensorize(const SizeType n_var, UVector &t) {
        ReadAndWrite<UVector> w(t);
        auto r = range(t);

        for (auto i = r.begin(); i < r.end(); i += n_var) {
            const auto value = t.get(i);

            for (SizeType k = 1; k < n_var; ++k) {
                t.set(i + k, value);
            }
        }
    }

}  // namespace utopia
