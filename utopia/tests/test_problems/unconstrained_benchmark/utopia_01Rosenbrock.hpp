#ifndef UTOPIA_ROSENBROCK_01
#define UTOPIA_ROSENBROCK_01

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_UnconstrainedTestFunction.hpp"


namespace utopia
{
    /**
     * @brief      Rosenbrock 2D banana function. \n 
     *             The floor of the valley follows approximately the parabola \f$ y = x^2 + 1/200 \f$.   
     *             The covariance matrix is not positive-definite. On the dashed line it is singular. 
     *             Stepping method tend to perform at least as well as gradient methods for this function.
     *
     */
    template<class Matrix, class Vector>
    class Rosenbrock01 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Rosenbrock01() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(this->dim());
            x_exact_ = values(this->dim(), 1.0);

            {
                const Write<Vector> write1(x_init_);
                x_init_.set(0, -1.2);
                x_init_.set(1, 1.0);        
            }

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == this->dim());

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result = 1 + 100.0 * pow(x * x - y , 2.0) + pow(x - 1 , 2.0);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            assert(point.size().get(0) == this->dim());
            result = zeros(this->dim());

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
            result.set(1, 200.0 * (y - x * x));
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
            assert(point.size().get(0) == this->dim());

            result = zeros(this->dim(), this->dim());

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar mixed = -400.0 * x;

            result.set(0, 0, 1200 * x * x - 400 * y + 2);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, 200.0);
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
            return 1; 
        }

        std::string name() const override
        {
            return "Rosenbrock"; 
        }

        SizeType dim() const override
        {
            return 2; 
        }


    private: 
        Vector x_init_; 
        Vector x_exact_; 

    };
}
#endif //UTOPIA_ROSENBROCK_01