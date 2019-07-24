#ifndef UTOPIA_GULF_11
#define UTOPIA_GULF_11

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Gulf11 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Gulf11()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(3);
            x_exact_ = zeros(3);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                x_init_.set(0, 5.0);
                x_init_.set(1, 2.5);
                x_init_.set(2, 0.15);

                // close one...
                // x_init_.set(0, 45.0);
                // x_init_.set(1, 22.0);
                // x_init_.set(2, 1.2);

                x_exact_.set(0, 50.0);
                x_exact_.set(1, 25.0);
                x_exact_.set(2, 1.5);
            }

        }

        std::string name() const override
        {
            return "Gulf reasearch and development";
        }

        SizeType dim() const override
        {
            return 3;
        }


        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 3);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            result = 0.0;
            for(SizeType i = 1; i <=99; i++)
            {
                Scalar arg = i*0.01;
                Scalar r = std::abs(std::pow(- 50.0* std::log(arg), 2.0/3.0) + 25.0 - y);
                Scalar t = std::exp(-std::pow(r, z)/x)-arg;

                result += t*t;
            }


            return true;
        }

        bool gradient(const Vector &point, Vector &g) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 3);
            g = zeros(3);

            const Read<Vector> read(point);
            const Write<Vector> write(g);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            Scalar g1 = 0.0;
            Scalar g2 = 0.0;
            Scalar g3 = 0.0;

            for(SizeType i =1; i <=99; i++)
            {
                Scalar arg = i*0.01;
                Scalar r = std::abs(std::pow(- 50.0* std::log(arg), 2.0/3.0) + 25.0 - y);

                Scalar t1 = std::pow(r, z)/x;
                Scalar t2 = std::exp(-t1);
                Scalar t = std::exp(-std::pow(r, z)/x)-arg;

                g1 += 2.0 * t * t1 * t2/x;
                g2 += 2.0 * t * t1 * t2* z / r;
                g3 -= 2.0 * t * t1 * t2 * std::log(r);
            }

            g.set(0, g1);
            g.set(1, g2);
            g.set(2, g3);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 3);
            result = zeros(3,3);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            Scalar term11 = 0.0;
            Scalar term22 = 0.0;
            Scalar term33 = 0.0;
            Scalar term21 = 0.0;
            Scalar term31 = 0.0;
            Scalar term32 = 0.0;

            for(SizeType i =1; i <=99; i++)
            {
                Scalar arg = i*0.01;
                Scalar r = std::pow( - 50.0 * std::log(arg), 2./3.) + 25.0 - y;
                Scalar t1 = std::pow(std::abs(r), z)/x;
                Scalar t2 = std::exp(-t1);
                Scalar t3 = t1 * t2 * ( t1 * t2 + (t1-1.0)*(t2-arg));
                Scalar t = t1 * t2 * (t2 - arg );
                Scalar logr = std::log(std::abs(r));

                term11 += 2.0 * t3 - (2.0 * t);
                term22 += 2.0 * (t+(z * t3))/r/r;
                term33 += 2.0 * t3 * logr * logr;
                term21 += 2.0 * t3 / r;
                term31 -= 2.0 * t3 * logr;
                term32 += 2.0 * (t-(z * t3 * logr))/r;
            }

            term11 = (term11/x)/x;
            term21 = term21*z/x;
            term22 = term22*z;
            term31 = term31/x;

            result.set(0, 0, term11);
            result.set(0, 1, term21);
            result.set(0, 2, term31);

            result.set(1, 0, term21);
            result.set(1, 1, term22);
            result.set(1, 2, term32);

            result.set(2, 0, term31);
            result.set(2, 1, term32);
            result.set(2, 2, term33);

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
            return 0.0;
        }

    private:
        Vector x_init_;
        Vector x_exact_;

    };
}

#endif //UTOPIA_GULF_11