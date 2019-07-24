#ifndef UTOPIA_BEALE_05
#define UTOPIA_BEALE_05

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Beale05 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Beale05()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = values(2, 1.0);
            x_exact_ = zeros(2);


            {
                const Write<Vector> write2(x_exact_);

                x_exact_.set(0, 3.0);
                x_exact_.set(1, 0.5);
            }

        }

        std::string name() const override
        {
            return "Beale";
        }

        SizeType dim() const override
        {
            return 2.0;
        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = 1.5 - x*(1.0 - y);
            Scalar b = 2.25 - x*(1.0 - y*y);
            Scalar c = 2.625 - x*(1.0 - y*y*y);

            result = a*a + b*b + c*c;
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }
            assert(point.size().get(0) == 2);
            result = zeros(2);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar y2 = std::pow(y, 2.0);
            const Scalar y3 = std::pow(y, 3.0);
            const Scalar y4 = std::pow(y, 4.0);
            const Scalar y5 = std::pow(y, 5.0);
            const Scalar y6 = std::pow(y, 6.0);

            Scalar a = 2.0*x * (y6 + y4 - 2.0*y3 - y2 - 2.0*y +3.0) + 5.25*y3 + 4.5 * y2 + 3.0*y - 12.75;
            Scalar b = x * (x * (6.0*y5 + 4.0*y3 - 6.0*y2 - 2.0*y - 2.0) + 15.75*y2 + 9*y +3);

            result.set(0, a);
            result.set(1, b);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }
            assert(point.size().get(0) == 2);

            result = zeros(2, 2);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar y2 = std::pow(y, 2.0);
            const Scalar y3 = std::pow(y, 3.0);
            const Scalar y4 = std::pow(y, 4.0);
            const Scalar y5 = std::pow(y, 5.0);
            const Scalar y6 = std::pow(y, 6.0);

            const Scalar mixed = 4.0 * x * (3.0 * y5 + 2.0 * y3 - 3.0 * y2 - y -1.0) + 15.75 * y2 + 9.0 * y + 3.0;
            const Scalar a = 2.0 * (y6 + y4 - 2.0*y3 - y2 - 2.0*y +3);
            const Scalar b = x * ( 2.0 * x * (15.0 * y4 + 6.0 * y2 - 6.0* y - 1.0) + 31.5 * y + 9);

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

#endif //UTOPIA_BEALE_05