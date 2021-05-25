#ifndef UTOPIA_PATCH_SMOOTHER_HPP
#define UTOPIA_PATCH_SMOOTHER_HPP

#include "utopia_PatchGatherer.hpp"

namespace utopia {

    template <class Matrix, class SerialMatrix = Matrix>
    class PatchSmoother {
    public:
        using Vector = typename Traits<Matrix>::Vector;
        using SerialVector = typename Traits<SerialMatrix>::Vector;
        using SizeType = typename Traits<Vector>::SizeType;
        using Scalar = typename Traits<Vector>::Scalar;

        bool apply(const Vector &b, Vector &x) {
            auto rr = b.range();

            Vector r, c;
            c.zeros(layout(x));

            r = *gatherer_.matrix() * x;
            r = b - r;

            Scalar norm_r = norm2(r);
            utopia::out() << "starting norm resdiual: " << norm_r << '\n';

            gatherer_.update(make_ref(r));

            for (int iter = 0; iter < 10000; ++iter) {
                for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                    gatherer_.local_to_patch_at_row(i);

                    if (!gatherer_.solve()) {
                        assert(false);
                        Utopia::Abort("Unable to solve patch!");
                    }

                    gatherer_.patch_to_local_at_row(i, c);
                }

                x += c;
                r = *gatherer_.matrix() * x;
                r = b - r;

                norm_r = norm2(r);

                utopia::out() << iter << ") norm resdiual: " << norm_r << '\n';
            }

            return true;
        }

        void update(const std::shared_ptr<const Matrix> &op) { gatherer_.update(op); }

        PatchSmoother() {
            gatherer_.set_patch_solver(std::make_shared<ConjugateGradient<SerialMatrix, SerialVector>>());
        }

        void set_patch_solver(const std::shared_ptr<LinearSolver<SerialMatrix, SerialVector>> &patch_solver) {
            gatherer_.set_patch_solver(patch_solver);
        }

    private:
        PatchGatherer<Matrix, SerialMatrix> gatherer_;
    };
}  // namespace utopia

#endif  // UTOPIA_PATCH_SMOOTHER_HPP
