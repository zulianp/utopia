#ifndef UTOPIA_SOLVER_TESTFUNCTIONSND_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONSND_HPP

#include <vector>
#include "utopia_Function.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class TestFunctionND_1 final : public Function<Matrix, Vector> {

        DEF_UTOPIA_SCALAR(Matrix)
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

    public:

        TestFunctionND_1(SizeType N = 10) : N(N), b(values(N, 3.0)), A(identity(N, N)), a(1.0) {

            static_assert(is_dense_or_polymorphic<Matrix>::value, "This function has a dense hessian do not use sparse implementations");
        }

        inline bool initialize_hessian(Matrix &H, Matrix &H_pre) const override
        {
            H = values(N, N, 0.);
            H_pre = values(N, N, 0.);
            return true;
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override {
            result =  0.5 * std::pow(dot(point, A * point) + 1, 2.0) - dot(b, point);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            if(empty(result) || size(point) != size(result)) {
                result = values(point.size(), 0.0);
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
        const Matrix A;
        const Scalar a;
    };

    template<class Matrix, class Vector>
    class SimpleQuadraticFunction : public Function<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        using SizeType = typename Traits<Vector>::SizeType;

        SimpleQuadraticFunction(const SizeType &n) : n_(n) {}

        virtual bool initialize_hessian(Matrix &H, Matrix &H_pre) const override
        {
            H = identity(n_, n_);
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
            result = identity(n, n);
            result *= 2;
            return true;
        }

        // SimpleQuadraticFunction() { }

        private:
            SizeType n_;
    };

    // Quadratic function class
    template<class Matrix, class Vector>
    class QuadraticFunctionBoundary : public Function<Matrix, Vector> {
    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        virtual bool value(const Vector &point, typename Vector::Scalar &result) const override {
            Scalar val = dot(point, A * point);
            Scalar val2 = dot(point, b);
            result = 0.5 * val - val2;
            return true;
        }

        virtual bool gradient(const Vector &point, Vector &result) const override {
            result = B*(A * point - b);
            //result = B * result;
            return true;
        }

        virtual bool hessian(const Vector &/*point*/, Matrix &result) const override {
            result = A;
            return true;
        }

        QuadraticFunctionBoundary(Vector b, Matrix H, Matrix B): b(b), A(H), B(B) { }
        private:
            const Vector b;
            Matrix A; /*!< Hessian */
            Matrix B; /*!< boundary operator */
    };
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Quadratic function class
    template<class Matrix, class Vector>
    class QuadraticFunctionConstrained : public ConstrainedTestFunction<Matrix, Vector> {
    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        QuadraticFunctionConstrained(Vector & b, Matrix & H, Matrix & B, Vector & ub):
                                                                        b_(b),
                                                                        A_(H),
                                                                        B_(B),
                                                                        ub_(ub)
        { 
            // to be figured out
            exact_sol_ = zeros(1); 
        }

        virtual bool value(const Vector &point, typename Vector::Scalar &result) const override {
            Scalar val = dot(point, A_ * point);
            Scalar val2 = dot(point, b_);
            result = 0.5 * val - val2;
            return true;
        }

        virtual bool gradient(const Vector &point, Vector &result) const override {
            result = B_*(A_ * point - b_);
            //result = B * result;
            return true;
        }

        virtual bool hessian(const Vector &/*point*/, Matrix &result) const override {
            result = A_;
            return true;
        }

        virtual bool upper_bound(Vector & ub) const override{
            ub = ub_;
            return true;
        }

        virtual bool lower_bound(Vector &/*lb*/) const override {
            assert("utopia::QuadraticFunctionConstrained:: This function doesnt have lower bound, just upper. \n");
            return false;
        }

        virtual bool has_upper_bound() const override
        {
            return true;
        }

        virtual bool has_lower_bound() const override
        {
            return false;
        }


        virtual Vector initial_guess() const override
        {   
            return (0 * b_);
        }
        
        virtual const Vector & exact_sol() const override
        {
            return exact_sol_; 
        }
        

        virtual Scalar min_function_value() const override
        {   
            // not known
            return 0.0; 
        }

        virtual std::string name() const override
        {
            return "Quadratic";
        }
        
        virtual SizeType dim() const override
        {
            return size(b_); 
        }

        virtual bool exact_sol_known() const override
        {
            return false;
        }

        virtual bool parallel() const override
        {
            return true;
        }

        private:
            const Vector b_;
            Matrix A_; /*!< Hessian */
            Matrix B_; /*!< boundary operator */
            const Vector ub_;
            Vector exact_sol_; 
    };


////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class Matrix, class Vector>
    class MildStiffExample : public virtual Function<Matrix, Vector> , public virtual LeastSquaresFunction<Matrix, Vector>
    {
        static_assert(!utopia::is_sparse<Matrix>::value || utopia::is_polymorhic<Matrix>::value, "utopia::MildStiffExample does not support sparse matrices as Hessian is dense matrix.");

    public:
        typedef UTOPIA_SCALAR(Vector)      Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)   SizeType;

        MildStiffExample(const SizeType & n): n_(n)
        {
            x_init_ = values(n_, 1.0);

            const SizeType n_local = local_size(x_init_);
            b_ = local_values(n_local, 1.0);
            Vector u = local_values(n_local, 1.0);

            Matrix U = outer(u, u);
            Scalar udot = 2./dot(u,u);
            Matrix I = local_identity(n_local, n_local);
            U = I - (udot * U);

            Matrix D = local_identity(n_local, n_local);

            {
                Write<Matrix> re(D);
                auto r = row_range(D);

                for(auto i=r.begin(); i != r.end(); ++i)
                    D.set(i,i, i+1);
            }

            // because some problem with petsc, when using UDU
            UDU_ = U * D;
            UDU_ *= U;

        }

        bool value(const Vector &x, Scalar &result) const override
        {
            assert(x.size() == n_);
            Vector g = 0*x;
            gradient(x, g);
            result = 0.5 * dot(g, g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(x.size() == n_);

            if(empty(g)){
                g = 0*x;
            }

            {
                Write<Vector> wg(g);
                Read<Vector> rx(x);
                auto r = range(g);

                for(auto i=r.begin(); i!=r.end(); ++i)
                    g.set(i, std::pow(x.get(i), 3.));
            }

            g = (UDU_* g) - b_;

            return true;
        }

        bool residual(const Vector &x, Vector &g) const override
        {
            return gradient(x, g);
        }

        bool jacobian(const Vector &x, Matrix &H) const override
        {
            return hessian(x, H);
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {

            Vector c = 0*x;

            {
                Write<Vector> wg(c);
                Read<Vector> rx(x);
                auto r = range(c);

                for(auto i=r.begin(); i!=r.end(); ++i)
                    c.set(i, std::pow(x.get(i), 2.));
            }

            Matrix C = diag(c);
            H = 3. * UDU_ * C;

            return true;
        }


        void get_initial_guess(Vector & x) const
        {
            x = x_init_;
        }

        private:
            const SizeType n_;
            Matrix UDU_;
            Vector b_;
            Vector x_init_;

    };


}


#endif //UTOPIA_SOLVER_TESTFUNCTIONSND_HPP
