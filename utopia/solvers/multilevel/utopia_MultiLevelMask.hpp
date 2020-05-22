#ifndef UTOPIA_MULTILEVEL_MASK_HPP
#define UTOPIA_MULTILEVEL_MASK_HPP

#include "utopia_Level.hpp"
#include "utopia_MatrixTransfer.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_Recorder.hpp"

#include <iostream>

namespace utopia {

    template <class Matrix, class Vector>
    class MultiLevelMask final {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::Transfer<Matrix, Vector> Transfer;
        using TransferPtr = std::shared_ptr<Transfer>;

    public:
        MultiLevelMask() = default;
        ~MultiLevelMask() = default;

        inline bool active() const { return active_; }

        void active(const bool val) { active_ = val; }

        void describe(std::ostream &os = std::cout) const {
            os << "active: " << active_ << "\n";
            os << "n_masks: " << masks_.size() << "\n";
        }

        void apply(const SizeType l, Vector &v) const {
            if (masks_.empty() || !active_) return;
            v = e_mul(masks_[l], v);
        }

        void generate_masks(const Matrix &A, const std::vector<TransferPtr> &transfers) {
            if (!active_) return;

            const auto L = transfers.size() + 1;

            masks_.resize(L);
            auto &mask = masks_[L - 1];

            generate_mask_from_matrix(A, mask, 0., 1.);
            // UTOPIA_RECORD_VALUE("r_mask", mask);

            for (SizeType l = L - 1; l > 0; --l) {
                auto &mask_l = masks_[l - 1];
                transfers[l - 1]->boolean_restrict_or(masks_[l], mask_l);

                // UTOPIA_RECORD_VALUE("r_mask", mask_l);
            }

            for (auto &m : masks_) {
                m.shift(-1);
                m = abs(m);
                // UTOPIA_RECORD_VALUE("mask", m);
            }
        }

    private:
        std::vector<Vector> masks_;
        bool active_{false};

        static void generate_mask_from_matrix(const Matrix &A,
                                              Vector &mask,
                                              const Scalar on_value,
                                              const Scalar off_value) {
            const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            // FIXME once atomic_store is avaialable
            mask.values(row_layout(A), off_value);

            // {
            //     auto mask_view = view_device(mask);

            //     A.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &value) {
            //         if (i == j) return;

            //         if (device::abs(value) > off_diag_tol) {
            //             mask_view.atomic_set(i, on_value);
            //         }
            //     });
            // }

            // mask.zeros(row_layout(A));

            {
                auto mask_view = view_device(mask);

                A.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &value) {
                    if (i == j) return;

                    if (device::abs(value) > off_diag_tol) {
                        // mask_view.atomic_add(i, 1.0);
                        mask_view.set(i, on_value);
                    }
                });
            }

            // mask.transform_values(UTOPIA_LAMBDA(const Scalar &v) {
            //     if (v >= 1.0) {
            //         return on_value;
            //     } else {
            //         return off_value;
            //     }
            // });

            // disp(mask);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MULTILEVEL_MASK_HPP
