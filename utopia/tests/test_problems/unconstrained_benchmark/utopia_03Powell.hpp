#ifndef UTOPIA_POWELL_03
#define UTOPIA_POWELL_03

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_UnconstrainedTestFunction.hpp"


namespace utopia
{
   
 template<class Matrix, class Vector>
    class Powell03 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        Powell03() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(2);
            x_exact_ = zeros(2);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                x_init_.set(0, 0.0);
                x_init_.set(1, 1.0);

                x_exact_.set(0, 1.09815933e-5);
                x_exact_.set(1, 9.106146738);                
            }

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = (10e4 * x * y -1.0); 
            Scalar b = std::exp(-x) + std::exp(-y) - 1.0001; 

            result = a*a + b*b; 
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            assert(point.size().get(0) == 2);
            result = zeros(2);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = 20000*(10000 *x * y -1.0); 
            Scalar b = std::exp(-x) + std::exp(-y) - 1.0001; 

            result.set(0, (y * a) - (2.0 * std::exp(-x) * b));
            result.set(1, (x * a) - (2.0 * std::exp(-y) * b));

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
            assert(point.size().get(0) == 2);

            result = zeros(2, 2);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            const Scalar mixed = 20000 * (20000*x * y - 1.0) + 2.*std::exp(-x -y); 
            const Scalar a = 2.*std::exp(-x-y) + 4.* std::exp(-2.0*x) - 2.0002*std::exp(-x) + 200000000*y*y; 
            const Scalar b = 200000000*x*x + 2.*std::exp(-x-y) + 4.* std::exp(-2.0*y) - 2.0002*std::exp(-y); 

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