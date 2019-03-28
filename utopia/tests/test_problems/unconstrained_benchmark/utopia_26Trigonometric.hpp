#ifndef UTOPIA_SOLVER_TRIGONOMETRIC_26
#define UTOPIA_SOLVER_TRIGONOMETRIC_26

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Trigonometric26 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Trigonometric26(const SizeType & n_loc): n_loc_(n_loc)
        {
            x_exact_ = local_values(n_loc_, 0.0);
            SizeType n_global = size(x_exact_).get(0);

            x_init_ = local_values(n_loc_, 1./Scalar(n_global));
            x_inc_ = local_values(n_loc_, 1.0);

            {
                const Write<Vector> write2(x_inc_);

                each_write(x_inc_, [](const SizeType i) -> double
                {
                    return i+1;
                }   );
            }
        }

        std::string name() const override
        {
            return "Trigonometric";
        }

        SizeType dim() const override
        {
            return n_loc_;
        }

        bool exact_sol_known() const override
        {
            return false;
        }


        bool value(const Vector &x, Scalar &result) const override
        {
            assert(local_size(x).get(0) == this->dim());

            SizeType n_global = size(x).get(0);

            Vector xcos = x;
            Vector xsin = x;

            {
                const Write<Vector> write1(xcos);
                const Write<Vector> write2(xsin);
                const Read<Vector> r1(x);

                each_write(xcos, [&x](const SizeType i) -> double
                {
                    return std::cos(x.get(i));
                }   );

                each_write(xsin, [&x](const SizeType i) -> double
                {
                    return std::sin(x.get(i));
                }   );
            }


            Scalar s1 = sum(xcos);
            Vector t = (n_global- s1) * local_values(local_size(x).get(0), 1.0);
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            result  = dot(t,t);

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(local_size(x).get(0) == this->dim());

            SizeType n_global = size(x).get(0);
            Vector xcos = x;
            Vector xsin = x;

            {
                const Write<Vector> write1(xcos);
                const Write<Vector> write2(xsin);
                const Read<Vector> r1(x);

                each_write(xcos, [&x](const SizeType i) -> double
                {
                    return std::cos(x.get(i));
                }   );

                each_write(xsin, [&x](const SizeType i) -> double
                {
                    return std::sin(x.get(i));
                }   );
            }

            Scalar s1 = sum(xcos);
            Vector t = (n_global- s1) * local_values(local_size(x).get(0), 1.0);
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            Scalar s2 = sum(t);
            g = e_mul(x_inc_, xsin) - xcos;
            g = e_mul(g, t);
            g = 2.0 * (g +(xsin*s2));

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(local_size(x).get(0) == this->dim());

            if(empty(H)){
                H = local_values(local_size(x).get(0), local_size(x).get(0), 0.0);
            }

            SizeType n_global = size(x).get(0);
            Vector xcos = x;
            Vector xsin = x;
            Vector ones = local_values(local_size(x).get(0), 1.0);

            {
                const Write<Vector> write1(xcos);
                const Write<Vector> write2(xsin);
                const Read<Vector> r1(x);

                each_write(xcos, [&x](const SizeType i) -> double
                {
                    return std::cos(x.get(i));
                }   );

                each_write(xsin, [&x](const SizeType i) -> double
                {
                    return std::sin(x.get(i));
                }   );
            }

            Scalar s1 = sum(xcos);
            Vector t = (n_global- s1) * ones;
            t = t + x_inc_ - xsin;
            t = t - e_mul(x_inc_, xcos);

            Scalar s2 = sum(t);

            {
                const Read<Vector> read(xsin);
                const Read<Vector> read2(xcos);
                const Read<Vector> read3(t);
                const Write<Matrix> write(H);

                each_write(H, [&xsin, &xcos, &t, s2, n_global](const SizeType i, const SizeType j) -> double
                {
                    Scalar val;

                    if(i==j)
                    {
                        val = ( ((i+1) * ( (i+1) + 2. )) + n_global ) * xsin.get(i) *xsin.get(i);
                        val += xcos.get(i) *( xcos.get(i) -  ( (2. * (i+1.)) + 2. ) * xsin.get(i));
                        val +=  t.get(i) * ( ((i+1.) * xcos.get(i)) + xsin.get(i));
                        val = 2.0 * ( val + (xcos.get(i) * s2));
                    }
                    else
                    {
                        Scalar th =  xcos.get(i);
                        Scalar xsj = xsin.get(i);
                        val = xsin.get(j) * ((( n_global + (i+1) + (j+1) ) * xsj) - th);
                        val -= (xsj * xcos.get(j));
                        val *= 2.0;
                    }

                    return val;

                }   );

            }

            return true;
        }

        Vector initial_guess() const override
        {
            return x_init_;
        }

        const Vector & exact_sol() const override
        {
            return x_exact_;
        }

        Scalar min_function_value() const override
        {
            return 0;
        }

    private:
        SizeType n_loc_;
        Vector x_init_;
        Vector x_exact_;
        Vector x_inc_;

    };

}

#endif //UTOPIA_SOLVER_TRIGONOMETRIC_26
