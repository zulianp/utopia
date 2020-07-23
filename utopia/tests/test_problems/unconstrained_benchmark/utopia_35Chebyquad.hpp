#ifndef UTOPIA_SOLVER_CHEBYQUAD_35
#define UTOPIA_SOLVER_CHEBYQUAD_35

#include <cassert>
#include <vector>
#include "utopia_Function.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class Chebyquad35 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Chebyquad35() { init(); }

        void init() {
            auto v_layout = serial_layout(dim());

            x_exact_.zeros(v_layout);
            x_init_.zeros(v_layout);
            SizeType n_global = 8;

            {
                // Device side writing
                auto x_view = view_device(x_init_);
                parallel_for(
                    range_device(x_init_),
                    UTOPIA_LAMBDA(const SizeType &i) { x_view.set(i, (i + 1) / Scalar(n_global + 1)); });

                // Host side writing
                Write<Vector> w(x_exact_);

                x_exact_.set(0, 0.043153);
                x_exact_.set(1, 0.193091);
                x_exact_.set(2, 0.266329);
                x_exact_.set(3, 0.500000);
                x_exact_.set(4, 0.500000);
                x_exact_.set(5, 0.733671);
                x_exact_.set(6, 0.806910);
                x_exact_.set(7, 0.956847);
            }
        }

        std::string name() const override { return "Chebyquad"; }

        SizeType dim() const override { return 8; }

        bool value(const Vector &x, Scalar &result) const override {
            if (x.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(size(x).get(0) == this->dim());

            Vector fvec;
            this->eval_polynomial(x, fvec);
            result = dot(fvec, fvec);

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override {
            if (x.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(size(x).get(0) == this->dim());

            Vector fvec;
            this->eval_polynomial(x, fvec);

            g.zeros(layout(x));
            std::vector<Scalar> g_help(this->dim());

            {
                const Read<Vector> read1(x);
                const Read<Vector> read2(fvec);

                for (auto j = 0; j < this->dim(); j++) {
                    Scalar t1 = 1.0;
                    Scalar t2 = (2.0 * x.get(j)) - 1.0;
                    Scalar t = 2.0 * t2;
                    Scalar s1 = 0.0;
                    Scalar s2 = 2.0;

                    for (auto i = 0; i < this->dim(); i++) {
                        g_help[j] += fvec.get(i) * s2;
                        Scalar th = (4.0 * t2) + (t * s2) - s1;
                        s1 = s2;
                        s2 = th;
                        th = (t * t2) - t1;
                        t1 = t2;
                        t2 = th;
                    }
                }
            }

            {
                const Write<Vector> write(g);

                for (auto i = 0; i < this->dim(); i++) {
                    g.set(i, g_help[i] / this->dim() * 2.0);
                }
            }

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override {
            if (x.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            auto n = this->dim();
            assert(size(x).get(0) == n);

            Vector fvec;
            this->eval_polynomial(x, fvec);

            H.dense(square_matrix_layout(layout(x)), 0.0);
            std::vector<std::vector<Scalar> > hess(n, std::vector<Scalar>(n));
            std::vector<Scalar> g(n);

            Scalar d1 = 1.0 / n;
            Scalar d2 = 2.0 * d1;

            {
                const Read<Vector> read1(fvec);
                const Read<Vector> read2(x);

                for (auto j = 1; j <= n; j++) {
                    hess[j - 1][j - 1] = 4.0 * d1;
                    Scalar t1 = 1.0;
                    Scalar t2 = (2.0 * x.get(j - 1)) - 1.0;
                    Scalar t = 2.0 * t2;
                    Scalar s1 = 0.0;
                    Scalar s2 = 2.0;
                    Scalar p1 = 0.0;
                    Scalar p2 = 0.0;
                    g[j - 1] = s2;

                    for (auto i = 2; i <= n; i++) {
                        Scalar th = (4.0 * t2) + (t * s2) - s1;
                        s1 = s2;
                        s2 = th;
                        th = (t * t2) - t1;
                        t1 = t2;
                        t2 = th;
                        th = (8.0 * s1) + (t * p2) - p1;
                        p1 = p2;
                        p2 = th;
                        g[i - 1] = s2;
                        hess[j - 1][j - 1] += (fvec.get(i - 1) * th) + (d1 * s2 * s2);
                    }

                    hess[j - 1][j - 1] *= d2;

                    for (auto k = 1; k <= j - 1; k++) {
                        hess[j - 1][k - 1] = 0.0;
                        Scalar tt1 = 1.0;
                        Scalar tt2 = (2.0 * x.get(k - 1)) - 1.0;
                        Scalar tt = 2.0 * tt2;
                        Scalar ss1 = 0.0;
                        Scalar ss2 = 2.0;

                        for (auto i = 1; i <= n; i++) {
                            hess[j - 1][k - 1] += ss2 * g[i - 1];
                            Scalar tth = (4.0 * tt2) + (tt * ss2) - ss1;
                            ss1 = ss2;
                            ss2 = tth;
                            tth = (tt * tt2) - tt1;
                            tt1 = tt2;
                            tt2 = tth;
                        }

                        hess[j - 1][k - 1] = d2 * d1 * hess[j - 1][k - 1];
                    }
                }
            }

            {
                const Write<Matrix> wr1(H);

                for (auto i = 0; i < this->dim(); i++) {
                    for (auto j = 0; j < this->dim(); j++) {
                        H.set(i, j, hess[i][j]);
                    }
                }
            }

            // this could be done way much nicer...
            Matrix D = diag(diag(H));
            H = H + transpose(H) - D;

            return true;
        }

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

        Scalar min_function_value() const override {
            return 3.51687e-3;  // if n = 10
        }

    private:
        void eval_polynomial(const Vector &x, Vector &fvec) const {
            if (empty(fvec)) {
                // fvec=local_zeros(local_size(x).get(0));
                fvec.zeros(layout(x));
            }

            const SizeType n_global = size(x).get(0);
            std::vector<Scalar> f_help(n_global);

            {
                const Read<Vector> read1(x);

                for (auto j = 0; j < n_global; j++) {
                    Scalar t1 = 1.0;
                    Scalar t2 = 2.0 * x.get(j) - 1.0;
                    Scalar t = 2.0 * t2;

                    for (auto i = 0; i < n_global; i++) {
                        f_help[i] = f_help[i] + t2;
                        Scalar th = t * t2 - t1;
                        t1 = t2;
                        t2 = th;
                    }
                }
            }

            for (auto i = 1; i <= n_global; i++) {
                f_help[i - 1] = f_help[i - 1] / n_global;

                if ((i % 2) == 0) {
                    f_help[i - 1] = f_help[i - 1] + (1.0 / ((i * i) - 1.0));
                }
            }

            {
                const Write<Vector> write(fvec);

                for (auto i = 0; i < n_global; i++) {
                    fvec.set(i, f_help[i]);
                }
            }
        }

    private:
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_CHEBYQUAD_35
