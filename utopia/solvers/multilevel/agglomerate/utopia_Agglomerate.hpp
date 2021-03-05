#ifndef UTOPIA_AGGLOMERATE_HPP
#define UTOPIA_AGGLOMERATE_HPP

#include "utopia_AlgebraicMultigrid.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class AdditiveCorrectionTransfer final : public Transfer<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using IndexArray = typename utopia::Traits<Vector>::IndexArray;

        bool interpolate(const Vector &x_coarse, Vector &x_fine) const override {
            auto vl = layout(x_coarse);
            assert(vl.local_rows() == n_coarse_local_);

            assert(!empty(x_fine));

            auto x_coarse_view = local_view_device(x_coarse);
            auto x_fine_view = local_view_device(x_fine);

            for (SizeType i = 0; i < n_coarse_local_; ++i) {
                x_fine_view.set(i, x_coarse_view.get(parent_[i]));
            }

            return false;
        }

        bool restrict(const Vector &x_fine, Vector &x_coarse) const override {
            auto vl = layout(x_fine);

            const SizeType n_fine = vl.rows();

            if (empty(x_coarse)) {
                x_coarse.zeros(layout(x_fine.comm(), n_coarse_local_, n_coarse_global_));
            } else {
                x_coarse.set(0.0);
            }

            auto x_fine_view = local_view_device(x_fine);
            auto x_coarse_view = local_view_device(x_coarse);
            for (SizeType i = 0; i < n_fine; ++i) {
                x_coarse_view.add(parent_[i], x_fine_view.get(i));
            }

            return true;
        }

        bool restrict(const Matrix & /*M_fine*/, Matrix & /*M_coarse*/) const override {
            assert(false);

            // IndexArray d_nnz(n_coarse_local_, 0);
            // // FIXME when doing a proper parallel version
            // auto &o_nnz = d_nnz;

            // PetscCrsView crs_fine(M_fine.raw_type());

            // std::vector<bool> touched(n_coarse_local_, false);

            // const SizeType n_fine = M_fine.local_rows();

            // for (SizeType i = 0; i < n_fine; ++i) {
            //     const auto &&row = crs_fine.row(i);
            //     const SizeType n_values = row.n_blocks();
            //     const SizeType parent_i = parent_[i];

            //     for (SizeType k = 0; k < n_values; ++k) {
            //         const SizeType j = row.colidx(k);
            //         const SizeType parent_j = parent_[j];

            //         if (!touched[parent_j]) {
            //             touched[parent_j] = true;
            //             ++d_nnz[parent_i];
            //         }
            //     }

            //     for (SizeType k = 0; k < n_values; ++k) {
            //         const SizeType j = row.colidx(k);
            //         const SizeType parent_j = parent_[j];
            //         touched[parent_j] = false;
            //     }
            // }

            // if (BlockSize > 1) {
            //     assert(false && "IMPLEMENT ME");
            //     return false;
            // } else {
            // M_coarse.sparse(layout(M_fine.comm(), n_coarse_local_, n_coarse_local_, n_coarse_global_,
            // n_coarse_global_),
            //                 d_nnz,
            //                 o_nnz);
            // }

            return false;
        }

        bool boolean_restrict_or(const Vector &, Vector &) override {
            assert(false && "IMPLENT ME");
            return false;
        }

        bool project_down(const Vector &, Vector &) const override {
            assert(false && "IMPLENT ME");
            return false;
        }

        bool project_down_positive_negative(const Vector &, const Vector &, Vector &) override {
            assert(false && "IMPLENT ME");
            return false;
        }

        void init_memory() override {}
        Scalar interpolation_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }
        Scalar projection_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }
        Scalar restriction_inf_norm() const override {
            assert(false && "IMPLENT ME");
            return -1.0;
        }

        void handle_equality_constraints(const Vector &) override { assert(false && "IMPLENT ME"); }

        IndexArray &parent() { return parent_; }

        void set_size(const SizeType n_coarse_local, SizeType n_coarse_global) {
            n_coarse_local_ = n_coarse_local;
            n_coarse_global_ = n_coarse_global;
        }

    private:
        SizeType n_coarse_local_{-1}, n_coarse_global_{-1};
        IndexArray parent_;
    };

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
        }

        std::shared_ptr<Transfer> create_transfer(const Matrix &in) override {
            UTOPIA_TRACE_REGION_BEGIN("Agglomerate::create_prolongator");

            using namespace utopia;

            auto prolongator = std::make_shared<Matrix>();

            auto rr = row_range(in);

            Matrix in_local;
            local_block_view(in, in_local);

            IndexArray parent(rr.extent(), -1);
            SizeType n_coarse_rows = 0;

            ScalarArray a_max(rr.extent(), 0);

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

                SizeType n_not_aggr = 0;
                for (SizeType i = 0; i < rr.extent(); ++i) {
                    if (parent[i] != -1) continue;

                    parent[i] = n_coarse_rows++;
                    ++n_not_aggr;
                }
            }

            auto pl = layout(in.comm(), in.local_rows(), n_coarse_rows, in.rows(), Traits::determine());
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
                    prolongator->set(rr.begin() + i, coarse_offset + parent[i], 1.0);
                }
            }

            UTOPIA_TRACE_REGION_END("Agglomerate::create_prolongator");

            return std::make_shared<IPRTransfer<Matrix, Vector>>(prolongator);
        }

        inline void verbose(const bool val) { verbose_ = val; }

    private:
        SizeType bmax_{3};
        // SizeType bmin_{5};
        Scalar weight_{1. / 3};
        bool verbose_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_AGGLOMERATE_HPP
