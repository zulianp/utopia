#ifndef UTOPIA_POWELL_03
#define UTOPIA_POWELL_03

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_Communicator.hpp"

#include <cassert>


namespace utopia
{

    template<class Matrix, class Vector>
    class Powell03 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        Powell03()
        {
            auto v_layout = serial_layout(2);

            x_init_.zeros(v_layout);
            x_exact_.zeros(v_layout);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                x_init_.set(0, 0.0);
                x_init_.set(1, 1.0);

                x_exact_.set(0, 1.09815933e-5);
                x_exact_.set(1, 9.106146738);
            }

        }

        std::string name() const override
        {
            return "Powell badly scaled";
        }

        SizeType dim() const override
        {
            return 2.0;
        }


        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
            if(point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = ((10000.0 * x * y) -1.0);
            Scalar b = std::exp(-x) + std::exp(-y) - 1.0001;

            result = a*a + b*b;
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override
        {
            if(point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);
            result.zeros(layout(point));

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar f1 = (10000.0 * x * y) - 1.0;
            Scalar df1dx1 = 10000.0 * y;
            Scalar df1dx2 = 10000.0 * x;

            Scalar f2 = std::exp(-x) + std::exp(-y) - 1.0001;
            Scalar df2dx1 = -std::exp(-x);
            Scalar df2dx2 = -std::exp(-y);

            Scalar a = (2.0 * f1 * df1dx1) + (2.0 * f2 * df2dx1);
            Scalar b = (2.0 * f1 * df1dx2) + (2.0 * f2 * df2dx2);

            result.set(0, a);
            result.set(1, b);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            if(point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 2);

            result.dense(square_matrix_layout(layout(point)), 0.0);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            const Scalar f1 = (10000.0 * x * y) - 1.0;
            const Scalar df1dx1 = 10000.0 * y;
            const Scalar df1dx2 = 10000.0 * x;
            const Scalar d2f1dx21 = 10000.0;

            const Scalar f2 = std::exp(-x) + std::exp(-y) - 1.0001;
            const Scalar df2dx1 = -std::exp(-x);
            const Scalar df2dx2 = -std::exp(-y);

            const Scalar d2f2dx11 = std::exp(-x);
            const Scalar d2f2dx22 = std::exp(-y);

            const Scalar mixed = (2.0 * f1 * d2f1dx21) + (2.0 * df1dx2 * df1dx1) + (2.0 * df2dx2 * df2dx1);
            const Scalar a = (2.0 * df1dx1 * df1dx1) + (2.0 * f2 * d2f2dx11) + (2.0 * df2dx1 * df2dx1);
            const Scalar b = (2.0 * df1dx2 * df1dx2) + (2.0 * f2 * d2f2dx22) + (2.0 * df2dx2 * df2dx2);

            result.set(0, 0, a);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, b);
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
        Vector x_init_;
        Vector x_exact_;
    };
}

#endif //UTOPIA_POWELL_03