#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Clonable.hpp"
#include "utopia_Trace.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Views.hpp"

#include <type_traits>
#include <vector>

namespace utopia {

    template <class Matrix>
    class ProjectedGaussSeidelSweep : public Clonable {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;
        using ArrayView = utopia::ArrayView<Scalar>;
        using ConstArrayView = utopia::ArrayView<const Scalar>;
        using VectorView = utopia::VectorView<ArrayView>;
        using ConstVectorView = utopia::VectorView<ConstArrayView>;

        virtual ~ProjectedGaussSeidelSweep() = default;
        virtual void init_from_local_matrix(const Matrix &mat) = 0;
        virtual void update_from_local_matrix(const Matrix &local_diag_block) = 0;
        virtual void apply(const SizeType &times) = 0;
        virtual void apply_unconstrained(const SizeType &times) = 0;
        ProjectedGaussSeidelSweep *clone() const override = 0;

        void set_residual_view(const ConstVectorView &r) {
            // FIXME
            r_.raw_type().set(&r[0], r.size());
        }

        void set_correction_view(VectorView &c) {
            // FIXME
            c_.raw_type().set(&c[0], c.size());
        }

        void set_correction_view(VectorView &&c) {
            // FIXME
            c_.raw_type().set(&c[0], c.size());
        }

        void set_bounds(const ConstVectorView &lb, const ConstVectorView &ub) {
            // FIXME
            lb_.raw_type().set(&lb[0], lb.size());
            ub_.raw_type().set(&ub[0], ub.size());
        }

        void symmetric(const bool value) { symmetric_ = value; }

        void l1(const bool l1) { l1_ = l1; }

        ProjectedGaussSeidelSweep() = default;

    protected:
        VectorView c_;
        ConstVectorView r_, lb_, ub_;

        bool symmetric_{true};
        bool l1_{false};
    };

    template <class Matrix>
    class ProjectedScalarGaussSeidelSweep final : public ProjectedGaussSeidelSweep<Matrix> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;
        using ArrayView = utopia::ArrayView<Scalar>;
        using ConstArrayView = utopia::ArrayView<const Scalar>;
        using VectorView = utopia::VectorView<ArrayView>;
        using ConstVectorView = utopia::VectorView<ConstArrayView>;

        ProjectedScalarGaussSeidelSweep *clone() const override { return new ProjectedScalarGaussSeidelSweep(); }

        void init_from_local_matrix(const Matrix &mat) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidelSweep::init_from_local_matrix");

            const SizeType n_rows = mat.rows();

            if (n_rows != SizeType(d_inv_.size())) {
                d_inv_.resize(n_rows);
                row_ptr_.resize(n_rows + 1);
            }

            std::fill(row_ptr_.begin(), row_ptr_.end(), SizeType(0));

            SizeType n_off_diag_entries = 0;

            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                if (i == j) {
                    d_inv_[i] = (device::abs(a_ij) > 0.0) ? (1 / a_ij) : Scalar(0.0);
                } else {
                    if (a_ij != 0.0) {
                        ++row_ptr_[i + 1];
                        ++n_off_diag_entries;
                    }
                }
            });

            if (this->l1_) {
                mat.read(
                    [&](const SizeType &i, const SizeType &, const Scalar &a_ij) { d_inv_[i] += device::abs(a_ij); });
            }

            for (SizeType i = 0; i < n_rows; ++i) {
                row_ptr_[i + 1] += row_ptr_[i];
            }

            assert(n_off_diag_entries == row_ptr_.back());

            if (n_off_diag_entries != SizeType(values_.size())) {
                values_.resize(n_off_diag_entries);
                col_idx_.resize(n_off_diag_entries);
            }

            n_off_diag_entries = 0;
            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                if (i != j && a_ij != 0.0) {
                    col_idx_[n_off_diag_entries] = j;
                    values_[n_off_diag_entries] = a_ij;
                    ++n_off_diag_entries;
                }
            });

            UTOPIA_TRACE_REGION_END("ProjectedGaussSeidelSweep::init_from_local_matrix");
        }

        void update_from_local_matrix(const Matrix &local_diag_block) override {
            init_from_local_matrix(local_diag_block);
        }

        void apply(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidelSweep::apply(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply();
            }

            UTOPIA_TRACE_REGION_END("ProjectedGaussSeidelSweep::apply(...)");
        }

        void apply_unconstrained(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidelSweep::apply_unconstrained(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply_unconstrained();
            }

            UTOPIA_TRACE_REGION_END("ProjectedGaussSeidelSweep::apply_unconstrained(...)");
        }

        ProjectedScalarGaussSeidelSweep() = default;

    private:
        void apply_unconstrained() {
            const SizeType n_rows = d_inv_.size();

            for (SizeType i = 0; i < n_rows; ++i) {
                Scalar val = this->r_[i];

                const SizeType row_end = row_ptr_[i + 1];
                for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                    val -= values_[j] * this->c_[col_idx_[j]];
                }

                this->c_[i] = d_inv_[i] * val;
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType i = n_rows - 1; i >= 0; --i) {
                    Scalar val = this->r_[i];

                    const SizeType row_end = row_ptr_[i + 1];
                    for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val -= values_[j] * this->c_[col_idx_[j]];
                    }

                    this->c_[i] = d_inv_[i] * val;
                }
            }
        }

        void apply() {
            const SizeType n_rows = d_inv_.size();

            for (SizeType i = 0; i < n_rows; ++i) {
                Scalar val = this->r_[i];

                const SizeType row_end = row_ptr_[i + 1];
                for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                    val -= values_[j] * this->c_[col_idx_[j]];
                }

                this->c_[i] = device::max(this->lb_[i], device::min(d_inv_[i] * val, this->ub_[i]));
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType i = n_rows - 1; i >= 0; --i) {
                    Scalar val = this->r_[i];

                    const SizeType row_end = row_ptr_[i + 1];
                    for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val -= values_[j] * this->c_[col_idx_[j]];
                    }

                    this->c_[i] = device::max(this->lb_[i], device::min(d_inv_[i] * val, this->ub_[i]));
                }
            }
        }

        std::vector<Scalar> values_;
        std::vector<Scalar> d_inv_;
        std::vector<SizeType> row_ptr_;
        std::vector<SizeType> col_idx_;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP
