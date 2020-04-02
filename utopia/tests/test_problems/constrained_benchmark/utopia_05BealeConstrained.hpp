#ifndef UTOPIA_BEALE_05_CONSTRAINED
#define UTOPIA_BEALE_05_CONSTRAINED

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Beale05Constrained final: public ConstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Beale05Constrained()
        {
            auto v_layout = serial_layout(dim());

            Vector ub, lb;
            ub.zeros(v_layout);
            lb.zeros(v_layout);

            {
                const Write<Vector> write1(ub);
                const Write<Vector> write2(lb);

                lb.set(0, 0.6);
                lb.set(1, 0.5);

                ub.set(0, 10);
                ub.set(1, 100);
            }


            this->set_box_constraints(make_box_constaints(std::make_shared<Vector>(lb), std::make_shared<Vector>(ub)));

        }

        std::string name() const override
        {
            return "Beale, bound constrained";
        }

        SizeType dim() const override
        {
            return unconstrained_.dim();
        }

        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            return unconstrained_.value(x, result);
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            return unconstrained_.gradient(x, g);
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            return unconstrained_.hessian(x, H);
        }

        Vector initial_guess() const override
        {
            return unconstrained_.initial_guess();
        }

        const Vector & exact_sol() const override
        {
            return unconstrained_.exact_sol();
        }

        Scalar min_function_value() const override
        {
            return 0.0;
        }

    private:
        Beale05<Matrix, Vector> unconstrained_;

    };
}

#endif //UTOPIA_BEALE_05_CONSTRAINED