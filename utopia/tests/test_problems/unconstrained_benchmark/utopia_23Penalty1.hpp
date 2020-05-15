#ifndef UTOPIA_SOLVER_PENALTY1_23
#define UTOPIA_SOLVER_PENALTY1_23

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class PenaltyI23 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        PenaltyI23(const SizeType &n_loc) : n_loc_(n_loc) { init(Comm::get_default(), n_loc); }

        PenaltyI23(const Comm &comm = Comm::get_default(), const SizeType &n_loc = 10) : n_loc_(n_loc) {
            init(comm, n_loc);
        }

        void init(const Comm &comm, const SizeType &n_loc) {
            x_init_.zeros(layout(comm, n_loc, Traits::determine()));
            x_exact_.values(layout(x_init_), 0.15812);  // depends on size.. this is valid for n=10

            {
                const Write<Vector> write1(x_init_);
                each_write(x_init_, [](const SizeType i) -> double { return i + 1; });
            }
        }

        std::string name() const override { return "Penalty I"; }

        SizeType dim() const override { return n_loc_; }

        bool parallel() const override { return true; }

        bool exact_sol_known() const override { return false; }

        bool value(const Vector &x, Scalar &result) const override {
            assert(local_size(x).get(0) == this->dim());

            Scalar alpha = 0.00001;
            Scalar t1 = -0.25 + dot(x, x);

            // Vector help = x - local_values(local_size(x).get(0), 1.0);
            Vector help = x;
            help.shift(-1.0);

            Scalar t2 = dot(help, help);

            result = (alpha * t2) + (t1 * t1);

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override {
            assert(local_size(x).get(0) == this->dim());

            Scalar alpha = 0.00001;
            Scalar t1 = -0.25 + dot(x, x);
            // Vector help = x - local_values(local_size(x).get(0), 1.0);
            Vector help = x;
            help.shift(-1.0);

            g = 2.0 * alpha * help;
            g += 4.0 * t1 * x;

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override {
            assert(local_size(x).get(0) == this->dim());

            H = outer(x, x);
            H *= 8.0;

            {
                const Read<Vector> read(x);
                const Write<Matrix> write(H);

                Scalar alpha = 0.00001;
                Scalar t1 = -0.25 + dot(x, x);
                Scalar d = (2.0 * alpha) + (4.0 * t1);

                auto r = row_range(H);
                for (auto i = r.begin(); i != r.end(); ++i) {
                    H.set(i, i, d + 8.0 * x.get(i) * x.get(i));
                }
            }

            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        Scalar min_function_value() const override {
            return 7.08765e-5;  // if n=10
        }

    private:
        SizeType n_loc_;
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_PENALTY1_23
