#ifndef UTOPIA_PATCH_SMOOTHER_HPP
#define UTOPIA_PATCH_SMOOTHER_HPP

#include "utopia_petsc_PatchGatherer.hpp"

namespace utopia {

    template <class Matrix, class SerialMatrix = Matrix>
    class PatchSmoother : public QPSolver<Matrix, typename Traits<Matrix>::Vector> {
    public:
        using Vector = typename Traits<Matrix>::Vector;
        using SerialVector = typename Traits<SerialMatrix>::Vector;
        using SizeType = typename Traits<Vector>::SizeType;
        using Scalar = typename Traits<Vector>::Scalar;
        using Super = utopia::QPSolver<Matrix, Vector>;

        inline PatchSmoother *clone() const override { return new PatchSmoother(*this); }

        bool apply(const Vector &b, Vector &x) override {
            auto rr = b.range();

            Vector r, c, u;
            c.zeros(layout(x));

            r = *gatherer_.matrix() * x;
            r = b - r;

            this->fill_empty_bounds(layout(b));

            Vector lb = *this->lower_bound();
            Vector ub = *this->upper_bound();

            ub -= x;
            lb -= x;

            gatherer_.update_rhs(make_ref(r));
            gatherer_.update_sol(make_ref(c));
            gatherer_.update_bounds(make_ref(lb), make_ref(ub));

            bool converged = false;

            if (this->verbose()) {
                this->init_solver("PatchSmoother", {" it. ", "|| u - u_old ||"});
            }

            for (int iter = 0; iter < this->max_it(); ++iter) {
                for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                    gatherer_.local_to_patch_at_row(i);

                    if (!gatherer_.solve()) {
                        assert(false);
                        Utopia::Abort("Unable to solve patch!");
                    }

                    gatherer_.patch_to_local_at_row(i);
                }

                x += c;
                r = *gatherer_.matrix() * x;
                r = b - r;

                ub -= c;
                lb -= c;

                Scalar diff = norm2(c);

                c.set(0.0);

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
            gatherer_.update(op);
        }

        PatchSmoother() {
            gatherer_.set_patch_solver(std::make_shared<ProjectedGaussSeidel<SerialMatrix, SerialVector>>());
        }

        void set_patch_solver(const std::shared_ptr<QPSolver<SerialMatrix, SerialVector>> &patch_solver) {
            gatherer_.set_patch_solver(patch_solver);
        }

    private:
        PatchGatherer<Matrix, SerialMatrix> gatherer_;
    };
}  // namespace utopia

#endif  // UTOPIA_PATCH_SMOOTHER_HPP
