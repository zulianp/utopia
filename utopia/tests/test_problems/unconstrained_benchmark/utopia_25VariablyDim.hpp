#ifndef UTOPIA_SOLVER_VARIABLY_DIMENSIONED_25
#define UTOPIA_SOLVER_VARIABLY_DIMENSIONED_25

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class VariablyDim25 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        VariablyDim25(const SizeType & n_loc): n_loc_(n_loc)
        {
            x_init_ = local_zeros(n_loc_);
            x_exact_ = local_values(n_loc_, 1.0);
            x_inc_ = local_values(n_loc_, 1.0);

            SizeType n_global = size(x_exact_).get(0);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_inc_);

                each_write(x_init_, [n_global](const SizeType i) -> double
                {
                    return (n_global - i - 1.)/Scalar(n_global);
                }   );

                each_write(x_inc_, [](const SizeType i) -> double
                {
                    return i+1;
                }   );
            }
        }

        std::string name() const override
        {
            return "Variably dimensioned";
        }

        SizeType dim() const override
        {
            return n_loc_;
        }

        bool parallel() const override
        {
            return true;
        }

        bool value(const Vector &x, Scalar &result) const override
        {
            assert(local_size(x).get(0) == this->dim());

            Vector help = x - local_values(local_size(x).get(0), 1.0);

            Scalar f1 = dot(x_inc_, help);
            Scalar f11 = f1*f1;
            Scalar f2 = dot(help, help);

            result = f11 * (1.0 + f11) + f2;

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(local_size(x).get(0) == this->dim());

            Vector help = x - local_values(local_size(x).get(0), 1.0);
            Scalar f1 = dot(x_inc_, help);

            g = ((2.0 * f1) + (4.0*f1*f1*f1)) * x_inc_;
            g  = g + (2.0 * help);

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(local_size(x).get(0) == this->dim());

            Vector help = x - local_values(local_size(x).get(0), 1.0);
            Scalar f1 = dot(x_inc_, help);

            H = outer(x_inc_, x_inc_);
            H *= 2.0 + (12.0 * f1 * f1);

            Vector d = diag(H);
            {
                const Read<Vector> read(d);
                const Write<Matrix> write(H);

                auto r = row_range(H);
                for(auto i = r.begin(); i != r.end(); ++i)
                {
                    H.set(i,i, d.get(i) + 2.0);
                }
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

#endif //UTOPIA_SOLVER_VARIABLY_DIMENSIONED_25
