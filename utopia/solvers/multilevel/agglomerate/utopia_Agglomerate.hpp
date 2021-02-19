#ifndef UTOPIA_AGGLOMERATE_HPP
#define UTOPIA_AGGLOMERATE_HPP

#include "utopia_AlgebraicMultigrid.hpp"

namespace utopia {

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class Agglomerate final : public MatrixAgglomerator<Matrix> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Comm = typename Traits::Communicator;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexArray = typename Traits::IndexArray;
        using Vector = typename Traits::Vector;

        Agglomerate *clone() const override { return new Agglomerate(); }

        void create_prolongator(const Matrix &in, Matrix &prolongator) override {
            using namespace utopia;

            Vector am;
            in.row_abs_max(am);
            auto rr = row_range(in);

            Matrix in_local;
            local_block_view(in, in_local);

            IndexArray count(rr.extent(), 0);
            IndexArray fine_to_coarse_idx(rr.extent(), 0);

            std::vector<bool> aggregated(rr.extent(), false);

            SizeType n_coarse_rows = 0;

            {
                auto am_view = const_local_view_device(am);

                SizeType last_row = -1;
                SizeType coarse_row = -1;

                in_local.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                    if (aggregated[i]) return;

                    if (last_row != i) {
                        if (last_row >= 0) aggregated[last_row] = true;

                        last_row = i;
                        ++coarse_row;
                    }

                    auto a_max = am_view.get(i);
                    auto &count_i = count[i];

                    if (i != j) {
                        if (weight * a_max < device::abs(a_ij) && count_i < max_agg && !aggregated[j]) {
                            aggregated[j] = true;
                            ++count_i;
                            fine_to_coarse_idx[j] = coarse_row;
                        }
                    } else {
                        fine_to_coarse_idx[i] = coarse_row;

                        ++count_i;
                        ++n_coarse_rows;
                    }
                });
            }

            auto pl = layout(in.comm(), in.local_rows(), n_coarse_rows, in.rows(), Traits::determine());
            prolongator.sparse(pl, 1, 1);

            {
                Write<Matrix> w(prolongator);

                auto n = rr.extent();
                auto coarse_offset = prolongator.col_range().begin();

                for (SizeType i = 0; i < n; ++i) {
                    prolongator.set(rr.begin() + i, coarse_offset + fine_to_coarse_idx[i], 1.0);
                }
            }
        }

        SizeType max_agg{3};
        Scalar weight{1. / 3};
    };

}  // namespace utopia

#endif  // UTOPIA_AGGLOMERATE_HPP
