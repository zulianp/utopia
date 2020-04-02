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
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        VariablyDim25(const SizeType &n_loc)
        : n_loc_(n_loc)
        {
            init(Comm::get_default(), n_loc);
        }

        VariablyDim25(const Comm &comm = Comm::get_default(), const SizeType &n_loc = 10): n_loc_(n_loc)
        {
            init(comm, n_loc);
        }

        void init(const Comm &comm, const SizeType &n_loc)
        {
            //determine global size once and reuse everywhere
            x_init_.zeros(layout(comm, n_loc, Traits::determine()));
            auto x_layout = layout(x_init_);

            x_exact_.values(x_layout, 1.0);
            x_inc_.values(x_layout, 1.0);

            help_ = make_unique<Vector>(x_layout, 0.0);
            ones_.values(x_layout, 1.0);

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

            *help_ = x - ones_;

            Scalar f1 = dot(x_inc_, *help_);
            Scalar f11 = f1*f1;
            Scalar f2 = dot(*help_, *help_);

            result = f11 * (1.0 + f11) + f2;

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(local_size(x).get(0) == this->dim());

            *help_ = x - ones_;
            Scalar f1 = dot(x_inc_, *help_);

            g = ((2.0 * f1) + (4.0*f1*f1*f1)) * x_inc_;
            g  = g + (2.0 * (*help_));

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(local_size(x).get(0) == this->dim());

            *help_ = x - ones_;
            Scalar f1 = dot(x_inc_, *help_);

            H = outer(x_inc_, x_inc_);
            H *= 2.0 + (12.0 * f1 * f1);

            *help_ = diag(H);
            {
                const Read<Vector> read(*help_);
                const Write<Matrix> write(H);

                auto r = row_range(H);
                for(auto i = r.begin(); i != r.end(); ++i)
                {
                    H.set(i, i, help_->get(i) + 2.0);
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
        Vector ones_;
        std::unique_ptr<Vector> help_;

    };

}

#endif //UTOPIA_SOLVER_VARIABLY_DIMENSIONED_25
