#ifndef UTOPIA_RAS_PATCH_SMOOTHER_HPP
#define UTOPIA_RAS_PATCH_SMOOTHER_HPP

#include "utopia_petsc_PatchGatherer.hpp"

namespace utopia {

    template <class Matrix, class SerialMatrix = Matrix>
    class RASPatchSmoother : public QPSolver<Matrix, typename Traits<Matrix>::Vector> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using SerialTraits = utopia::Traits<SerialMatrix>;

        using Vector = typename Traits::Vector;
        using SerialVector = typename SerialTraits::Vector;
        using SizeType = typename Traits::SizeType;
        using Scalar = typename Traits::Scalar;
        using IndexArray = typename Traits::IndexArray;

        using Super = utopia::QPSolver<Matrix, Vector>;

        inline RASPatchSmoother *clone() const override { return new RASPatchSmoother(*this); }

        void read(Input &in) override {
            Super::read(in);
            in.get("overlap", overlap_);
            in.get("dumping", dumping_);
        }

        bool apply(const Vector &b, Vector &x) override {
            this->fill_empty_bounds(layout(b));

            Vector r, c, overlapping_r, overlapping_c, overlapping_ub, overlapping_lb;

            r = *this->get_operator() * x;
            r = b - r;

            Vector lb = *this->lower_bound();
            Vector ub = *this->upper_bound();

            ub -= x;
            lb -= x;

            overlapping_r = (*row_permutation_) * r;
            overlapping_lb = (*row_permutation_) * lb;
            overlapping_ub = (*row_permutation_) * ub;
            overlapping_c.zeros(layout(overlapping_r));

            gatherer_.update_rhs(make_ref(overlapping_r));
            gatherer_.update_sol(make_ref(overlapping_c));
            gatherer_.update_bounds(make_ref(overlapping_lb), make_ref(overlapping_ub));

            bool converged = false;

            if (this->verbose()) {
                this->init_solver("RASPatchSmoother", {" it. ", "|| u - u_old ||"});
            }

            auto rr = range(overlapping_c);

            for (int iter = 0; iter < this->max_it(); ++iter) {
                for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                    gatherer_.local_to_patch_at_row(i);

                    if (!gatherer_.solve()) {
                        assert(false);
                        Utopia::Abort("Unable to solve patch!");
                    }

                    gatherer_.patch_to_local_at_row(i);
                }

                c = transpose(*row_permutation_) * overlapping_c;
                c = e_mul(c, *weights_);
                c *= dumping_;

                x += c;
                r = *this->get_operator() * x;
                r = b - r;

                ub -= c;
                lb -= c;

                overlapping_r = (*row_permutation_) * r;
                overlapping_lb = (*row_permutation_) * lb;
                overlapping_ub = (*row_permutation_) * ub;

                Scalar diff = norm2(c);

                overlapping_c.set(0.0);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(iter, {diff});
                }

                if ((converged = this->check_convergence(iter, 1, 1, diff))) {
                    // disp(x);
                    break;
                }
            }

            return converged;
        }

        static void convert_to_sequential_matrix(const Matrix &parallel_matrix, Matrix &seq_matrix) {
            auto n_row_local = parallel_matrix.local_rows();
            auto n_cols = parallel_matrix.cols();

            IndexArray d_nnz(n_row_local, 0);
            IndexArray o_nnz(n_row_local, 0);

            auto rr = row_range(parallel_matrix);

            parallel_matrix.read(
                [&](const SizeType &i, const SizeType &, const Scalar &) { d_nnz[i - rr.begin()] += 1; });

            seq_matrix.sparse(serial_layout(n_row_local, n_cols), d_nnz, o_nnz);

            Write<Matrix> w(seq_matrix);

            parallel_matrix.read(
                [&](const SizeType &i, const SizeType &j, const Scalar &v) { seq_matrix.set(i - rr.begin(), j, v); });
        }

        static void convert_to_parallel_matrix(const Matrix &seq_matrix, Matrix &parallel_matrix) {
            auto n_local = seq_matrix.local_rows();

            IndexArray d_nnz(n_local, 0);
            IndexArray o_nnz(n_local, 0);

            seq_matrix.read([&](const SizeType &i, const SizeType &, const Scalar &) { d_nnz[i] += 1; });

            auto vl = layout(parallel_matrix.comm(), n_local, Traits::determine());
            parallel_matrix.sparse(square_matrix_layout(vl), d_nnz, o_nnz);

            auto rr = row_range(parallel_matrix);

            Write<Matrix> w(parallel_matrix);

            seq_matrix.read([&](const SizeType &i, const SizeType &j, const Scalar &v) {
                parallel_matrix.set(i + rr.begin(), j + rr.begin(), v);
            });
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);

            IS rows{nullptr}, cols{nullptr};
            MatGetOwnershipIS(op->raw_type(), &rows, &cols);
            MatIncreaseOverlap(op->raw_type(), 1, &rows, overlap_);

            auto local_rows = op->local_rows();
            auto global_rows = op->rows();

            PetscInt new_size;
            ISGetLocalSize(rows, &new_size);

            overlapping_matrix_ = std::make_shared<Matrix>(op->comm());
            row_permutation_ = std::make_shared<Matrix>();

            row_permutation_->sparse(layout(op->comm(), new_size, local_rows, Traits::determine(), global_rows), 1, 1);

            auto ml = layout(*row_permutation_);

            {
                Write<Matrix> w_row(*row_permutation_);
                // Write<Matrix> w_col(*col_permutation_);

                const PetscInt *index_map{nullptr};

                ISGetIndices(rows, &index_map);

                auto rr = row_range(*row_permutation_);

                for (PetscInt i = 0; i < new_size; ++i) {
                    row_permutation_->set(i + rr.begin(), index_map[i], 1.0);
                }

                ISRestoreIndices(rows, &index_map);
            }

            Vector ones(row_layout(*row_permutation_), 1.0);

            weights_ = std::make_shared<Vector>(row_layout(*op));
            *weights_ = transpose(*row_permutation_) * ones;
            (*weights_) = 1. / (*weights_);

            {
                Matrix temp = (*row_permutation_) * (*op);  // * (*col_permutation_);
                Matrix temp_seq;
                convert_to_sequential_matrix(temp, temp_seq);

                Matrix perm_seq;
                convert_to_sequential_matrix(*row_permutation_, perm_seq);

                Matrix overlapping_matrix_seq = temp_seq * transpose(perm_seq);
                convert_to_parallel_matrix(overlapping_matrix_seq, *overlapping_matrix_);
            }

            op->comm().barrier();

            // rename("op", const_cast<Matrix &>(*op));
            // write("load_op.m", *op);

            // rename("om", *overlapping_matrix_);
            // write("load_om.m", *overlapping_matrix_);

            // rename("perm", *row_permutation_);
            // write("load_perm.m", *row_permutation_);

            gatherer_.update(overlapping_matrix_);
        }

        RASPatchSmoother() {
            gatherer_.set_patch_solver(std::make_shared<ProjectedGaussSeidel<SerialMatrix, SerialVector>>());
        }

        void set_patch_solver(const std::shared_ptr<QPSolver<SerialMatrix, SerialVector>> &patch_solver) {
            gatherer_.set_patch_solver(patch_solver);
        }

        ~RASPatchSmoother() {}

        // void init_memory(const Layout &layout) {
        //     UTOPIA_TRACE_REGION_BEGIN("RASPatchSmoother::init_memory(...)");
        //     Super::init_memory(layout);

        //     UTOPIA_TRACE_REGION_END("RASPatchSmoother::init_memory(...)");
        // }

    private:
        PatchGatherer<Matrix, SerialMatrix> gatherer_;
        std::shared_ptr<Matrix> row_permutation_;
        std::shared_ptr<Vector> weights_;
        std::shared_ptr<Matrix> overlapping_matrix_;
        int overlap_{1};
        Scalar dumping_{2. / 3.};
    };
}  // namespace utopia

#endif  // UTOPIA_RAS_PATCH_SMOOTHER_HPP
