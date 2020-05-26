#ifndef UTOPIA_MILD_STIFF_EXAMPLE_HPP
#define UTOPIA_MILD_STIFF_EXAMPLE_HPP

#include <vector>
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class MildStiffExample : public virtual Function<Matrix, Vector>,
                             public virtual LeastSquaresFunction<Matrix, Vector> {
        static_assert(!utopia::is_sparse<Matrix>::value || utopia::is_polymorhic<Matrix>::value,
                      "utopia::MildStiffExample does not support sparse matrices as Hessian is dense matrix.");

    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        MildStiffExample(const SizeType &n) : n_(n) {
            auto x_layout = layout(Comm::get_default(), Traits::decide(), n);
            auto mat_layout = square_matrix_layout(x_layout);

            x_init_.values(x_layout, 1.0);

            b_.values(x_layout, 1.0);
            Vector u(x_layout, 1.0);

            Matrix U = outer(u, u);
            Scalar udot = 2. / dot(u, u);
            Matrix I;
            I.identity(mat_layout, 1.0);
            U = I - (udot * U);

            Matrix D;
            D.identity(mat_layout, 1.0);

            {
                Write<Matrix> re(D);
                auto r = row_range(D);

                for (auto i = r.begin(); i != r.end(); ++i) D.set(i, i, i + 1);
            }

            // because some problem with petsc, when using UDU
            UDU_ = U * D;
            UDU_ *= U;
        }

        bool value(const Vector &x, Scalar &result) const override {
            assert(x.size() == n_);
            Vector g = 0 * x;
            gradient(x, g);
            result = 0.5 * dot(g, g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override {
            assert(x.size() == n_);

            if (empty(g)) {
                g = 0 * x;
            }

            {
                Write<Vector> wg(g);
                Read<Vector> rx(x);
                auto r = range(g);

                for (auto i = r.begin(); i != r.end(); ++i) g.set(i, std::pow(x.get(i), 3.));
            }

            g = (UDU_ * g) - b_;

            return true;
        }

        bool residual(const Vector &x, Vector &g) const override { return gradient(x, g); }

        bool jacobian(const Vector &x, Matrix &H) const override { return hessian(x, H); }

        bool hessian(const Vector &x, Matrix &H) const override {
            Vector c = 0 * x;

            {
                Write<Vector> wg(c);
                Read<Vector> rx(x);
                auto r = range(c);

                for (auto i = r.begin(); i != r.end(); ++i) c.set(i, std::pow(x.get(i), 2.));
            }

            Matrix C = diag(c);
            H = 3. * UDU_ * C;

            return true;
        }

        void get_initial_guess(Vector &x) const { x = x_init_; }

    private:
        const SizeType n_;
        Matrix UDU_;
        Vector b_;
        Vector x_init_;
    };

}  // namespace utopia

#endif  // UTOPIA_MILD_STIFF_EXAMPLE_HPP
