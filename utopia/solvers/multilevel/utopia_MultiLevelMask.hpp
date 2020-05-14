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
        MultiLevelMask()

            = default;

        ~MultiLevelMask() = default;

        inline bool active() const { return active_; }

        void active(const bool val) { active_ = val; }

        // static void fix_semidefinite_operator(Matrix &A)
        // {

        // 	Vector d;

        // 	Size s = local_size(A);
        // 	d = local_values(s.get(0), 1.);

        // 	{
        // 		Write<Vector> w_d(d);

        // 		each_read(A,[&d](const SizeType i, const SizeType, const double) {
        // 			d.set(i, 0.);
        // 		});
        // 	}

        // 	A += Matrix(diag(d));
        // }

        // inline void set_fix_semidefinite_operators(const bool val)
        // {
        // 	fix_semidefinite_operators_ = val;
        // }

        void describe(std::ostream &os = std::cout) const {
            os << "active: " << active_ << "\n";
            os << "n_masks: " << masks_.size() << "\n";
        }

        void apply(const SizeType l, Vector &v) const {
            if (masks_.empty() || !active_) return;
            v = e_mul(masks_[l], v);
        }

        void generate_masks(const Matrix &A, const std::vector<TransferPtr> &transfers) {
            // const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;
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
                // m = local_values(local_size(m).get(0), 1.) - m;
                m.shift(-1);
                m = abs(m);
                // UTOPIA_RECORD_VALUE("mask", m);
            }
        }

    private:
        std::vector<Vector> masks_;
        // bool fix_semidefinite_operators_;
        bool active_{false};

        static void generate_mask_from_matrix(const Matrix &A,
                                              Vector &mask,
                                              const Scalar on_value,
                                              const Scalar off_value) {
            const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            // auto ls = local_size(A);
            mask.values(row_layout(A), off_value);

            {
                Write<Vector> w_(mask);

                each_read(A, [&](const SizeType i, const SizeType j, const Scalar value) {
                    if (i == j) return;

                    if (std::abs(value) > off_diag_tol) {
                        mask.set(i, on_value);
                    }
                });
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MULTILEVEL_MASK_HPP
