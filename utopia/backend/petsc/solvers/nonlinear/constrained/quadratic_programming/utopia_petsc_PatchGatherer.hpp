#ifndef UTOPIA_PATCH_GATHERER_HPP
#define UTOPIA_PATCH_GATHERER_HPP

#include "utopia_Layout.hpp"
#include "utopia_RowView.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class SerialMatrix = Matrix>
    class PatchGatherer : public Configurable {
    public:
        using Vector = typename Traits<Matrix>::Vector;
        using SerialVector = typename Traits<SerialMatrix>::Vector;
        using SizeType = typename Traits<Vector>::SizeType;
        using IndexArray = typename Traits<SerialVector>::IndexArray;

        class Patch {
        public:
            IndexArray col_idx;

            SerialMatrix patch_matrix_;
            SerialVector patch_rhs_;
            SerialVector patch_sol_;
            SerialVector patch_lb_;
            SerialVector patch_ub_;

            void patch_to_local_at_row(const Matrix &matrix, const SizeType row, Vector &x) {
                RowView<const Matrix> row_view(matrix, row);
                SizeType n_values = row_view.n_values();
                auto rr = row_range(matrix);

                Write<Vector> write_sol(x);
                Read<SerialVector> read_sol(patch_sol_);

                for (SizeType k = 0, patch_idx = 0; k < n_values; ++k) {
                    SizeType col = row_view.col(k);

                    if (rr.inside(col)) {
                        x.add(col, patch_sol_.get(patch_idx));
                        patch_idx++;
                    }
                }
            }

            void local_to_patch_at_row(const Vector &rhs, const Vector &lb, const Vector &ub, const Vector &x) {
                SizeType n_values = col_idx.size();
                SizeType n_patch = 0;

                auto rr = range(rhs);

                {
                    Read<Vector> read_rhs(rhs);
                    Read<Vector> read_x(x), read_lb(lb), read_ub(ub);

                    Write<SerialVector> write_rhs(patch_rhs_);
                    Write<SerialVector> write_x(patch_sol_), write_lb(patch_lb_), write_ub(patch_ub_);

                    for (SizeType k = 0, patch_row = 0; k < n_values; ++k) {
                        const SizeType col = col_idx[k];

                        if (rr.inside(col)) {
                            patch_rhs_.set(patch_row, rhs.get(col));
                            patch_sol_.set(patch_row, x.get(col));
                            patch_lb_.set(patch_row, lb.get(col));
                            patch_ub_.set(patch_row, ub.get(col));
                            ++patch_row;
                        }
                    }
                }

                patch_rhs_ -= patch_matrix_ * patch_sol_;
                patch_lb_ -= patch_sol_;
                patch_ub_ -= patch_sol_;

                patch_sol_.set(0.0);
            }

            void local_to_patch_matrix_at_row(IndexArray &matrix_to_patch, const Matrix &matrix, const SizeType r) {
                SizeType n_values = -1;
                SizeType n_patch = 0;
                auto rr = row_range(matrix);

                {
                    RowView<const Matrix> row(matrix, r);
                    n_values = row.n_values();

                    if (SizeType(col_idx.size()) < n_values) {
                        col_idx.resize(n_values);
                    }

                    if (matrix_to_patch.empty() || SizeType(matrix_to_patch.size()) < rr.extent()) {
                        matrix_to_patch.resize(rr.extent());
                        std::fill(matrix_to_patch.begin(), matrix_to_patch.end(), -1);
                    }

                    for (SizeType k = 0; k < n_values; ++k) {
                        SizeType col = row.col(k);

                        col_idx[k] = col;

                        if (rr.inside(col)) {
                            n_patch++;
                        }
                    }
                }

                auto local_layout = serial_layout(n_patch);
                auto matrix_layout = square_matrix_layout(local_layout);

                patch_rhs_.clear();
                patch_sol_.clear();

                patch_rhs_.zeros(local_layout);
                patch_sol_.zeros(local_layout);

                patch_lb_.zeros(local_layout);
                patch_ub_.zeros(local_layout);

                IndexArray d_nnz;
                IndexArray o_nnz;
                d_nnz.resize(n_patch);
                o_nnz.resize(n_patch);

                for (SizeType k = 0, patch_idx = 0; k < n_values; ++k) {
                    SizeType col = col_idx[k];

                    if (rr.inside(col)) {
                        assert(patch_idx < SizeType(d_nnz.size()));

                        RowView<const Matrix> other_row(matrix, col);

                        o_nnz[patch_idx] = 0;
                        // Assumes diagonal is there
                        const SizeType local_idx = col - rr.begin();

                        assert(local_idx >= 0);
                        assert(local_idx < rr.extent());

                        matrix_to_patch[local_idx] = patch_idx;

                        d_nnz[patch_idx++] = other_row.n_values();
                    }
                }

                assert(matrix_to_patch[r - rr.begin()] != -1);

                patch_matrix_.clear();
                patch_matrix_.sparse(matrix_layout, d_nnz, o_nnz);

                {
                    Write<SerialMatrix> write_matrix(patch_matrix_);

                    for (SizeType k = 0; k < n_values; ++k) {
                        const SizeType col = col_idx[k];

                        if (rr.inside(col)) {
                            RowView<const Matrix> other_row(matrix, col);
                            const SizeType other_n_values = other_row.n_values();
                            const SizeType patch_row = matrix_to_patch[col - rr.begin()];
                            assert(patch_row != -1);

                            for (SizeType other_k = 0; other_k < other_n_values; ++other_k) {
                                const SizeType other_col = other_row.col(other_k);

                                if (rr.inside(other_col)) {
                                    const SizeType patch_col = matrix_to_patch[other_col - rr.begin()];

                                    if (patch_col == -1) {
                                        continue;
                                    }

                                    const auto val = other_row.get(other_k);

                                    patch_matrix_.set(patch_row, patch_col, val);
                                }
                            }
                        }
                    }
                }

                for (auto &c : col_idx) {
                    if (rr.inside(c)) {
                        matrix_to_patch[c - rr.begin()] = -1;
                    }
                }
            }
        };

        void read(Input &in) override { in.get("store_patches", store_patches_); }

        void update(const std::shared_ptr<const Matrix> &matrix) {
            matrix_ = matrix;

            if (store_patches_) {
                SizeType n = matrix->local_rows();
                patches_.resize(n);

                for (SizeType i = 0; i < n; ++i) {
                    patches_[i].local_to_patch_matrix_at_row(matrix_to_patch, *matrix, i);
                }
            } else {
                // patches_.resize(1);
            }
        }

        void update_rhs(const std::shared_ptr<const Vector> &rhs) { rhs_ = rhs; }
        void update_sol(const std::shared_ptr<Vector> &x) { x_ = x; }

        void update_bounds(const std::shared_ptr<Vector> &lb, const std::shared_ptr<Vector> &ub) {
            assert(lb);
            assert(ub);

            lb_ = lb;
            ub_ = ub;
        }

        /// @brief Overwrites the values connect to the row
        void patch_to_local_at_row(const SizeType row) {
            current_patch_ = row;

            if (store_patches_) {
                patches_[row].patch_to_local_at_row(*matrix_, row, *x_);
            } else {
                RowView<const Matrix> row_view(*matrix_, row);
                SizeType n_values = row_view.n_values();
                auto rr = row_range(*matrix_);

                Write<Vector> write_sol(*x_);
                Read<SerialVector> read_sol(patch_sol_);

                for (SizeType k = 0, patch_idx = 0; k < n_values; ++k) {
                    SizeType col = row_view.col(k);

                    if (rr.inside(col)) {
                        x_->add(col, patch_sol_.get(patch_idx));
                        patch_idx++;
                    }
                }
            }
        }

        void local_to_patch_at_row(const SizeType r) {
            current_patch_ = r;

            if (store_patches_) {
                patches_[current_patch_].local_to_patch_at_row(*rhs_, *lb_, *ub_, *x_);
            } else {
                SizeType n_values = -1;
                SizeType n_patch = 0;
                auto rr = row_range(*matrix_);

                {
                    RowView<const Matrix> row(*matrix_, r);
                    n_values = row.n_values();

                    if (SizeType(col_idx_.size()) < n_values) {
                        col_idx_.resize(n_values);
                    }

                    if (matrix_to_patch.empty() || SizeType(matrix_to_patch.size()) < rr.extent()) {
                        matrix_to_patch.resize(rr.extent());
                        std::fill(matrix_to_patch.begin(), matrix_to_patch.end(), -1);
                    }

                    for (SizeType k = 0; k < n_values; ++k) {
                        SizeType col = row.col(k);

                        col_idx_[k] = col;

                        if (rr.inside(col)) {
                            n_patch++;
                        }
                    }
                }

                auto local_layout = serial_layout(n_patch);
                auto matrix_layout = square_matrix_layout(local_layout);

                patch_rhs_.clear();
                patch_sol_.clear();

                patch_rhs_.zeros(local_layout);
                patch_sol_.zeros(local_layout);

                patch_lb_.zeros(local_layout);
                patch_ub_.zeros(local_layout);

                d_nnz.resize(n_patch);
                o_nnz.resize(n_patch);

                // #ifndef NDEBUG
                // matrix_to_patch[r - rr.begin()] = -1;
                // #endif  // NDEBUG

                for (SizeType k = 0, patch_idx = 0; k < n_values; ++k) {
                    SizeType col = col_idx_[k];

                    if (rr.inside(col)) {
                        assert(patch_idx < SizeType(d_nnz.size()));

                        RowView<const Matrix> other_row(*matrix_, col);

                        o_nnz[patch_idx] = 0;
                        // Assumes diagonal is there
                        const SizeType local_idx = col - rr.begin();

                        assert(local_idx >= 0);
                        assert(local_idx < rr.extent());

                        matrix_to_patch[local_idx] = patch_idx;

                        d_nnz[patch_idx++] = other_row.n_values();
                    }
                }

                assert(matrix_to_patch[r - rr.begin()] != -1);

                patch_matrix_.clear();
                patch_matrix_.sparse(matrix_layout, d_nnz, o_nnz);

                {
                    Read<Vector> read_rhs(*rhs_);
                    Read<Vector> read_x(*x_), read_lb(*lb_), read_ub(*ub_);

                    Write<SerialMatrix> write_matrix(patch_matrix_);
                    Write<SerialVector> write_rhs(patch_rhs_);
                    Write<SerialVector> write_x(patch_sol_), write_lb(patch_lb_), write_ub(patch_ub_);

                    for (SizeType k = 0; k < n_values; ++k) {
                        const SizeType col = col_idx_[k];

                        if (rr.inside(col)) {
                            RowView<const Matrix> other_row(*matrix_, col);
                            const SizeType other_n_values = other_row.n_values();
                            const SizeType patch_row = matrix_to_patch[col - rr.begin()];
                            assert(patch_row != -1);

                            for (SizeType other_k = 0; other_k < other_n_values; ++other_k) {
                                const SizeType other_col = other_row.col(other_k);

                                if (rr.inside(other_col)) {
                                    const SizeType patch_col = matrix_to_patch[other_col - rr.begin()];

                                    if (patch_col == -1) {
                                        continue;
                                    }

                                    const auto val = other_row.get(other_k);

                                    patch_matrix_.set(patch_row, patch_col, val);
                                }
                            }

                            patch_rhs_.set(patch_row, rhs_->get(col));
                            patch_sol_.set(patch_row, x_->get(col));
                            patch_lb_.set(patch_row, lb_->get(col));
                            patch_ub_.set(patch_row, ub_->get(col));
                        }
                    }
                }

                patch_rhs_ -= patch_matrix_ * patch_sol_;
                patch_lb_ -= patch_sol_;
                patch_ub_ -= patch_sol_;

                patch_sol_.set(0.0);

                patch_solver_->set_box_constraints(
                    BoxConstraints<SerialVector>(make_ref(patch_lb_), make_ref(patch_ub_)));

                // write("patch_" + std::to_string(r) + ".m", patch_matrix_);

                for (auto &c : col_idx_) {
                    if (rr.inside(c)) {
                        matrix_to_patch[c - rr.begin()] = -1;
                    }
                }
            }

            // if (n_patch > max_patch_size_) {
            //     std::cout << "n_patch: " << n_patch << std::endl;
            //     std::cout << "nnz " << patch_matrix_.nnz() << "/" << (patch_matrix_.rows() * patch_matrix_.cols())
            //               << std::endl;
            // }

            // max_patch_size_ = std::max(max_patch_size_, n_patch);
        }

        bool solve() {
            assert(patch_solver_);

            if (!patch_solver_) return false;

            if (store_patches_) {
                patch_solver_->set_box_constraints(BoxConstraints<SerialVector>(
                    make_ref(patches_[current_patch_].patch_lb_), make_ref(patches_[current_patch_].patch_ub_)));

                patches_[current_patch_].patch_sol_.set(0.0);
                return patch_solver_->solve(patches_[current_patch_].patch_matrix_,
                                            patches_[current_patch_].patch_rhs_,
                                            patches_[current_patch_].patch_sol_);
            } else {
                patch_sol_.set(0.0);
                return patch_solver_->solve(patch_matrix_, patch_rhs_, patch_sol_);
            }
        }

        void set_patch_solver(const std::shared_ptr<QPSolver<SerialMatrix, SerialVector>> &patch_solver) {
            patch_solver_ = patch_solver;
        }

        inline std::shared_ptr<const Matrix> matrix() { return matrix_; }

        PatchGatherer(const PatchGatherer &other) : patch_solver_(other.patch_solver_->clone()) {}
        PatchGatherer() {}

        ~PatchGatherer() {}

    private:
        std::shared_ptr<const Matrix> matrix_;
        std::shared_ptr<const Vector> rhs_;
        std::shared_ptr<Vector> x_;

        std::shared_ptr<Vector> lb_;
        std::shared_ptr<Vector> ub_;

        IndexArray d_nnz;
        IndexArray o_nnz;
        IndexArray matrix_to_patch;

        SerialMatrix patch_matrix_;
        SerialVector patch_rhs_;
        SerialVector patch_sol_;
        SerialVector patch_lb_;
        SerialVector patch_ub_;

        std::shared_ptr<QPSolver<SerialMatrix, SerialVector>> patch_solver_;

        IndexArray col_idx_;

        SizeType max_patch_size_{0};
        SizeType current_patch_{0};

        bool store_patches_{true};
        std::vector<Patch> patches_;
    };
}  // namespace utopia

#endif  // UTOPIA_PATCH_GATHERER_HPP
