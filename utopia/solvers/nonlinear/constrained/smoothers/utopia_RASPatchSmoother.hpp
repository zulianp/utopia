#ifndef UTOPIA_RAS_PATCH_SMOOTHER_HPP
#define UTOPIA_RAS_PATCH_SMOOTHER_HPP

#include "utopia_PatchGatherer.hpp"

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

        using Super = utopia::QPSolver<Matrix, Vector>;

        inline RASPatchSmoother *clone() const override { return new RASPatchSmoother(*this); }

        void read(Input &in) override {
            Super::read(in);
            in.get("overlap", overlap_);
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

            overlapping_r = (*permutation_) * r;
            overlapping_lb = (*permutation_) * lb;
            overlapping_ub = (*permutation_) * ub;
            overlapping_c.zeros(layout(overlapping_r));

            gatherer_.update_rhs(make_ref(overlapping_r));
            gatherer_.update_sol(make_ref(overlapping_c));
            gatherer_.update_bounds(make_ref(overlapping_lb), make_ref(overlapping_ub));

            bool converged = false;

            if (this->verbose()) {
                this->init_solver("utopia::RASPatchSmoother", {" it. ", "|| u - u_old ||"});
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

                c = transpose(*permutation_) * overlapping_c;
                c = e_mul(c, *weights_);

                x += c;
                r = *this->get_operator() * x;
                r = b - r;

                ub -= c;
                lb -= c;

                overlapping_r = (*permutation_) * r;

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

        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);

            IS rows{nullptr}, cols{nullptr};
            MatGetOwnershipIS(op->raw_type(), &rows, &cols);
            MatIncreaseOverlap(op->raw_type(), 1, &rows, overlap_);

            auto local_rows = op->local_rows();
            auto global_rows = op->rows();

            PetscInt new_size;
            ISGetLocalSize(rows, &new_size);

            overlapping_matrix_ = std::make_shared<Matrix>();
            permutation_ = std::make_shared<Matrix>();

            permutation_->sparse(layout(op->comm(), new_size, local_rows, Traits::determine(), global_rows), 1, 1);

            {
                Write<Matrix> w(*permutation_);

                const PetscInt *index_map{nullptr};

                ISGetIndices(rows, &index_map);

                auto rr = row_range(*permutation_);

                for (PetscInt i = 0; i < new_size; ++i) {
                    permutation_->set(i + rr.begin(), index_map[i], 1.0);
                }

                ISRestoreIndices(rows, &index_map);
            }

            Vector ones(row_layout(*permutation_), 1.0);

            weights_ = std::make_shared<Vector>(row_layout(*op));
            *weights_ = transpose(*permutation_) * ones;
            (*weights_) = 1. / (*weights_);

            (*overlapping_matrix_) = (*permutation_) * (*op) * transpose(*permutation_);

            // disp(*overlapping_matrix_);
            // exit(0);

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
        std::shared_ptr<Matrix> permutation_;
        std::shared_ptr<Vector> weights_;
        std::shared_ptr<Matrix> overlapping_matrix_;
        int overlap_{1};
    };
}  // namespace utopia

#endif  // UTOPIA_RAS_PATCH_SMOOTHER_HPP
