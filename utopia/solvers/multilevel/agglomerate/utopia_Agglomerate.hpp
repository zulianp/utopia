#ifndef UTOPIA_AGGLOMERATE_HPP
#define UTOPIA_AGGLOMERATE_HPP

#include "utopia_AdditiveCorrectionTransfer.hpp"
#include "utopia_IPRTruncatedTransfer.hpp"
#include "utopia_MatrixAgglomerator.hpp"
#include "utopia_RowView.hpp"

namespace utopia {

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class Agglomerate final : public MatrixAgglomerator<Matrix> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Comm = typename Traits::Communicator;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexArray = typename Traits::IndexArray;
        using ScalarArray = typename Traits::ScalarArray;
        using Vector = typename Traits::Vector;
        using Transfer = utopia::Transfer<Matrix, Vector>;

        Agglomerate *clone() const override { return new Agglomerate(*this); }

        void read(Input &in) override {
            in.get("bmax", bmax_);
            in.get("weight", weight_);
            in.get("verbose", verbose_);
            in.get("use_additive_correction", use_additive_correction_);
        }

        std::shared_ptr<Transfer> create_transfer(const Matrix &in) override {
            return std::make_shared<IPRTransfer<Matrix, Vector>>(create_transfer_matrix(in));
        }

        std::shared_ptr<Transfer> create_truncated_transfer(const Matrix &in) override {
            return std::make_shared<IPRTruncatedTransfer<Matrix, Vector>>(create_transfer_matrix(in));
        }

        inline void verbose(const bool val) { verbose_ = val; }
        inline void use_additive_correction(const bool val) { use_additive_correction_ = val; }

    private:
        SizeType bmax_{3};
        // SizeType bmin_{5};
        Scalar weight_{1. / 3};
        bool verbose_{false};
        bool use_additive_correction_{false};

        std::shared_ptr<Matrix> create_transfer_matrix(const Matrix &in) {
            UTOPIA_TRACE_REGION_BEGIN("Agglomerate::create_transfer_matrix");

            using namespace utopia;

            auto rr = row_range(in);

            Matrix in_local;
            local_block_view(in, in_local);

            IndexArray parent(rr.extent(), -1);
            SizeType n_coarse_rows = 0;

            ScalarArray a_max(rr.extent(), 0);

            if (!use_additive_correction_) {
                for (SizeType i = 0; i < rr.extent(); ++i) {
                    RowView<const Matrix> row_view(in_local, i);
                    if (row_view.n_values() == 1) {
                        parent[i] = -2;
                    }
                }
            }

            // Global a_max_i
            // in.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
            //     if (i != j) {
            //         auto local_i = i - rr.begin();
            //         a_max[local_i] = device::max(a_max[local_i], device::abs(a_ij));
            //     }
            // });

            // Local a_max_i
            in_local.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                if (i != j) {
                    a_max[i] = device::max(a_max[i], device::abs(a_ij));
                }
            });

            {
                for (SizeType i = 0; i < rr.extent(); ++i) {
                    if (parent[i] != -1) continue;

                    parent[i] = n_coarse_rows;

                    RowView<const Matrix> row_view(in_local, i);
                    SizeType n_values = row_view.n_values();

                    Scalar a_max_i = a_max[i];

                    SizeType count = 0;
                    for (SizeType k = 0; k < n_values; ++k) {
                        const SizeType j = row_view.col(k);
                        const Scalar a_ij = row_view.get(k);

                        if (weight_ * a_max_i < device::abs(a_ij) && parent[j] == -1) {
                            parent[j] = n_coarse_rows;
                            count++;
                        }

                        if (count >= bmax_) {
                            break;
                        }
                    }

                    if (count == 0) {
                        // This node is not a cluster center
                        parent[i] = -1;

                        Scalar max_aij = 0;
                        SizeType arg_max_j = -1;
                        for (SizeType k = 0; k < n_values; ++k) {
                            const SizeType j = row_view.col(k);
                            const Scalar a_ij = std::abs(row_view.get(k));

                            if (a_ij > max_aij) {
                                arg_max_j = j;
                                max_aij = a_ij;
                            }
                        }

                        if (arg_max_j != -1) {
                            parent[i] = parent[arg_max_j];
                        }

                    } else {
                        ++n_coarse_rows;
                    }
                }

                for (SizeType i = 0; i < rr.extent(); ++i) {
                    if (parent[i] != -1) continue;

                    parent[i] = n_coarse_rows++;
                }
            }

            std::shared_ptr<Matrix> ret;

            // if (use_additive_correction_) {
            //     auto act = std::make_shared<AdditiveCorrectionTransfer<Matrix, Vector>>();
            //     act->parent() = std::move(parent);

            //     SizeType n_global_coarse_row = in.comm().sum(n_coarse_rows);
            //     act->set_size(n_coarse_rows, n_global_coarse_row);
            //     ret = act;
            // } else {
            auto pl = layout(in.comm(), in.local_rows(), n_coarse_rows, in.rows(), Traits::determine());

            auto prolongator = std::make_shared<Matrix>();
            prolongator->sparse(pl, 1, 1);

            if (verbose_) {
                in.comm().synched_print(
                    std::to_string(prolongator->rows()) + " -> " + std::to_string(prolongator->cols()) +
                    ", coarsening factor: " + std::to_string(prolongator->rows() / float(prolongator->cols())) + '\n');
            }

            {
                Write<Matrix> w(*prolongator);

                auto n = rr.extent();
                auto coarse_offset = prolongator->col_range().begin();

                for (SizeType i = 0; i < n; ++i) {
                    if (parent[i] >= 0) {
                        prolongator->set(rr.begin() + i, coarse_offset + parent[i], (a_max[i] != 0) ? 1 : 0);
                    }
                }
            }

            ret = prolongator;
            // }

            UTOPIA_TRACE_REGION_END("Agglomerate::create_transfer_matrix");
            return ret;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_AGGLOMERATE_HPP
