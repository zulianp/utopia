//
// Created by Alessandro Rigazzi on 28/05/15.
//


#ifndef UTOPIA_SOLVER_TESTFUNCTIONSND_HPP
#define UTOPIA_SOLVER_TESTFUNCTIONSND_HPP

#include <vector>
#include "utopia_Function.hpp"
#include "utopia_Core.hpp"


namespace utopia {

    template<class Matrix, class Vector>
    class TestFunctionND_1 final : public Function<Matrix, Vector> {

        DEF_UTOPIA_SCALAR(Matrix)
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

    public:

        TestFunctionND_1(SizeType N = 10) : N(N), b(values(N, 3.0)), A(identity(N, N)), a(1.0) {

            static_assert(Matrix::FILL_TYPE == FillType::DENSE, "This function has a dense hessian do not use sparse implementations");
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
            if(empty(result) || size(point).get(0) != size(result).get(0)) {
                result = values(point.size().get(0), 0.0);
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
            const SizeType n = size(point).get(0);

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
            const auto n = point.size().get(0);
            result = identity(n, n);
            result *= 2;
            return true;
        }

        SimpleQuadraticFunction() { }
    };

    // Quadratic function class
    template<class Matrix, class Vector>
    class QuadraticFunctionBoundary : public Function<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

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

        virtual bool hessian(const Vector &point, Matrix &result) const override {
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
    class QuadraticFunctionConstrained : public FunctionBoxConstrained<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        QuadraticFunctionConstrained(Vector & b, Matrix & H, Matrix & B, Vector & ub):
                                                                        b_(b),
                                                                        A_(H),
                                                                        B_(B),
                                                                        ub_(ub)
        { }

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

        virtual bool hessian(const Vector &point, Matrix &result) const override {
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

        private:
            const Vector b_;
            Matrix A_; /*!< Hessian */
            Matrix B_; /*!< boundary operator */
            const Vector ub_;
    };


////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Matrix, class Vector>
    class RosenbrockGeneric : public Function<Matrix, Vector>
    {
    public:
        typedef typename utopia::Traits<Vector>::SizeType SizeType;
        DEF_UTOPIA_SCALAR(Matrix)

        RosenbrockGeneric() { }

        bool update(const Vector &point) override
        {
            init_perm(point);
            x_perm_ = perm_ * point;

            {
                Read<Vector> read_x(x_perm_);
                Range r_x = range(x_perm_);
                xm1 = x_perm_.get(r_x.begin());
                xp1 = x_perm_.get(r_x.begin() + 1);
            }

            // Read<Vector> read_p(point);
            // Range r = r_p(point);
            //
            // // computing value
            // value_ = 0;
            // for (SizeType i = r.begin(); i < r.end(); i++) {
            //     Scalar xi = point.get(i);
            //     Scalar xnext = point.get(i+1);
            //     sum += 100.0 * pow(xnext - xi * xi, 2.0) + pow(xi - 1, 2.0);
            // }

            return true;
        }

        bool value(const Vector &point, Scalar &result) const override
        {
            const Read<Vector> read(point);
            const Range r = range(point);

            Scalar sum = 0;
            SizeType endm1 = r.end() - 1;
            for (SizeType i = r.begin(); i < endm1; i++) {
                Scalar xi = point.get(i);
                Scalar xnext = point.get(i+1);
                sum += 100.0 * pow(xnext - xi * xi, 2.0) + pow(xi - 1, 2.0);
            }
            if (r.end() != point.size().get(0))
                sum += 100.0 * pow(xp1 - point.get(endm1) * point.get(endm1), 2.0)
                    + pow(point.get(endm1) - 1, 2.0);
            result = sum;
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override
        {
            Read<Vector> read(point);

            SizeType d = point.size().get(0);
            result = zeros(d);

            // for parallel access to result
            Range r_result = range(result);
            if (r_result.empty())
                return false;
            SizeType rr_begin = r_result.begin();
            SizeType rr_end = r_result.end();

            // for parallel access to point
            Range r_point = range(point);
            SizeType rp_begin = r_point.begin();
            SizeType rp_end = r_point.end();

            // is this true?
            // assert(rr_begin == rp_begin && rr_end == rp_end);

            // must be put after range otherwise it will crash (segfault)
            Write<Vector> write(result);

            if (rr_begin == 0) {
                Scalar x0 = point.get(0);
                Scalar x1 = (1 == rp_end ? xp1 : point.get(1));
                result.set(0, - 400.0 * (x1 - x0 * x0) * x0 + 2.0 * (x0 - 1));
                rr_begin++;
            }

            if (rr_end == d) {
                Scalar dm2 = (d-1 == rp_begin ? xm1 : point.get(d-2));
                result.set(d-1, 200.0 * (point.get(d-1) - dm2 * dm2));
                rr_end--;
            }

            for (SizeType i = rr_begin; i < rr_end; i++) {
                Scalar xi = point.get(i);
                Scalar xprev = (i == rp_begin ? xm1 : point.get(i-1));
                Scalar xnext = (i+1 == rp_end ? xp1 : point.get(i+1));

                result.set(i,
                    200.0 * (xi - xprev * xprev)
                    - 400.0 * (xnext - xi * xi) * xi + 2.0 * (xi - 1)
                );
            }

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            Read<Vector> read(point);

            SizeType d = point.size().get(0);
            result = zeros(d, d);

            // for parallel access to Matrix
            Range r_result = row_range(result);
            if (r_result.empty())
                return false;
            SizeType r_begin = r_result.begin();
            SizeType r_end = r_result.end();

            // for parallel access to point
            Range r_point = range(point);
            SizeType rp_begin = r_point.begin();
            SizeType rp_end = r_point.end();

            Write<Matrix> write(result);

            // first row
            if (r_begin == 0 && r_begin != r_end) {
                Scalar x0 = point.get(0);
                Scalar x1 = (1 == rp_end ? xp1 : point.get(1));
                result.set(0, 0,
                    - 400.0 * ((x1 - x0 * x0) - 2.0 * x0 * x0) + 2.0
                );
                result.set(0, 1, - 400.0 * x0);
                r_begin++;
            }

            // last row
            if (r_end == d && r_begin != r_end) {
                result.set(d-1, d-2, - 400.0 * (d-1 == rp_begin ? xm1 : point.get(d-2)));
                result.set(d-1, d-1, 200.0);
                r_end--;
            }

            // all other rows
            for (SizeType i = r_begin; i < r_end; i++) {
                Scalar xi = point.get(i);
                Scalar xprev = (i == rp_begin ? xm1 : point.get(i-1));
                Scalar xnext = (i+1 == rp_end ? xp1 : point.get(i+1));

                result.set(i, i-1, - 400.0 * xprev);
                result.set(i, i,
                    200.0 - 400.0 * ((xnext - xi * xi) - 2.0 * xi * xi) + 2.0
                );
                result.set(i, i+1, - 400.0 * xi);
            }

            return true;
        }

        void clear()
        {

        }

    private:
        Matrix perm_;
        Vector x_perm_;
        Scalar xm1, xp1;

        Scalar value_;

        void init_perm(const Vector &x)
        {
            // if(!empty(perm_)) {
            //     return;
            // }

            auto r = range(x);
            long n = local_size(x).get(0);
            long N = size(x).get(0);

            if (is_sparse<Matrix>::value) {
                perm_ = local_sparse(2, n, 2);
            } else {
                perm_ = local_zeros({2, n});
            }

            auto r_perm = row_range(perm_);
            const long i1 = r_perm.begin();
            const long i2 = r_perm.begin() + 1;

            {
                Write<Matrix> w_p(perm_);

                if(r.begin() > 0) {
                    perm_.set(i1, r.begin() - 1, 1);
                }

                if(r.end() != N) {
                    perm_.set(i2, r.end(), 1);
                }
            }
            // disp(perm_);

        }
    };


    template<class Matrix, class Vector>
    class MildStiffExample : public Function<Matrix, Vector> 
    {
        static_assert(!utopia::is_sparse<Matrix>::value, "utopia::MildStiffExample does not support sparse matrices as Hessian is dense matrix.");

    public:
        typedef UTOPIA_SCALAR(Vector)      Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)   SizeType;

        MildStiffExample(const SizeType & n): n_(n)
        {
            x_init_ = values(n_, 0.0);    

            const SizeType n_local = local_size(x_init_).get(0); 
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
            assert(x.size().get(0) == 2);
            Vector g = values(2, 0.0); 
            gradient(x, g); 
            result = 0.5 * norm2(g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override 
        {
            g = local_values(local_size(x).get(0), 0.0);

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

        bool hessian(const Vector &x, Matrix &H) const override 
        {

            Vector c = local_values(local_size(x).get(0), 0.0);

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
