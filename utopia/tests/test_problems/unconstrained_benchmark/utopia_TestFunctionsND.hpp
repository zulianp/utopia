#ifndef UTOPIA_SOLVER_TESTFUNCTIONSND_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONSND_HPP

#include <vector>
#include "utopia_Function.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class TestFunctionND_1 final : public Function<Matrix, Vector> {
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

    public:

        TestFunctionND_1(SizeType N)
        : N(N),
          b(layout(Comm::get_default(), Traits::decide(), N), 3.0),
          A(layout(Comm::get_default(), Traits::decide(), Traits::decide(), N, N)),
          a(1.0)
          {
            A.identity();
            help_ = make_unique<Vector>(layout(b), 0.0);

            static_assert(is_dense_or_polymorphic<Matrix>::value, "This function has a dense hessian do not use sparse implementations");
        }

        TestFunctionND_1(const Comm &comm = Comm(), SizeType N = 10)
        : N(N),
          b(layout(comm, Traits::decide(), N), 3.0),
          A(layout(comm, Traits::decide(), Traits::decide(), N, N)),
          a(1.0)
          {
            A.identity();
            help_ = make_unique<Vector>(layout(b), 0.0);

            static_assert(is_dense_or_polymorphic<Matrix>::value, "This function has a dense hessian do not use sparse implementations");
        }


        inline bool initialize_hessian(Matrix &H, Matrix &H_pre) const override
        {
            H.dense(layout(A), 0.0);
            H_pre.dense(layout(A), 0.0);
            return true;
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            *help_ = A * point;
            result =  0.5 * std::pow(dot(point, *help_) + 1, 2.0) - dot(b, point);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if(empty(result) || size(point) != size(result)) {
                result.zeros(layout(point));
            }

            const Scalar s = dot(point, A * point);
            const Range r = range(point);

            {
                const Read<Vector> r_p(point);
                const Read<Vector> r_b(b);

                const Write<Vector> write(result);

                for (auto i = r.begin(); i != r.end(); ++i) {
                    result.set(i, 2.0 * a * point.get(i) * (s + 1.0) - b.get(i));
                }
            }

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            Matrix temp = outer(point, point);

            if(empty(result) || size(result) != size(temp)) {
                result = temp;
            }

            const Scalar s = dot(point, A * point);
            const Range rr = row_range(temp);
            const SizeType n = size(point);

            assert(rr.begin() == range(point).begin());
            assert(rr.end() == range(point).end());

            const Scalar identityScaleFactor = 2.0 * s + 2.0;

            {
                const Read<Matrix>  read(temp);
                const Write<Matrix> write(result);

                for (auto i = rr.begin(); i != rr.end(); ++i) {
                    for (SizeType j = 0; j != n; ++j) {
                        result.set(i, j, 4.0 * temp.get(i, j) + (SizeType(i) == j) * (identityScaleFactor));
                    }
                }
            }

            return true;
        }

        bool has_preconditioner() const override
        {
            return true;
        }

        bool hessian(const Vector &point, Matrix &H, Matrix &precond) const override
        {
            if(!hessian(point, H)) return false;

            if(empty(precond)) {
                precond = diag(diag(H));
            } else {
                Vector d = diag(H);
                precond *= 0.;

                Write<Matrix> w_p(precond);
                Read<Vector>  r_d(d);

                auto r = row_range(precond);

                for(auto i = r.begin(); i != r.end(); ++i) {
                    precond.set(i, i, d.get(i));
                }

            }

            return true;
        }

    private:
        const SizeType N;
        const Vector b;
        Matrix A;
        const Scalar a;
        std::unique_ptr<Vector>  help_;
    };

    template<class Matrix, class Vector>
    class SimpleQuadraticFunction : public Function<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        using SizeType = typename Traits<Vector>::SizeType;

        SimpleQuadraticFunction(const SizeType &n) : n_(n) {}

        virtual bool initialize_hessian(Matrix &H, Matrix &H_pre) const override
        {
            H.identity(serial_layout(n_, n_), 2.0);
            H_pre = H;
            return true;
        }

        virtual bool value(const Vector &point, Scalar &result) const override {
            const Scalar val = norm2(point);
            result = val * val;
            return true;
        }

        virtual bool gradient(const Vector &point, Vector &result) const override {
            result = 2.0 * point;
            return true;
        }

        virtual bool hessian(const Vector &point, Matrix &result) const override {
            const auto n = point.size();
            result.identity(square_matrix_layout(layout(point)), 2.0);
            return true;
        }

        // SimpleQuadraticFunction() { }

        private:
            SizeType n_;
    };

}


#endif //UTOPIA_SOLVER_TESTFUNCTIONSND_HPP
