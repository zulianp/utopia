#ifndef UTOPIA_GAUSSIAN_09
#define UTOPIA_GAUSSIAN_09

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Gaussian09 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Gaussian09()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(3);
            x_exact_ = zeros(3);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                x_init_.set(0, 0.4);
                x_init_.set(1, 1.0);
                x_init_.set(2, 0.0);

                x_exact_.set(0, 0.398956);
                x_exact_.set(1, 1.00002);
                x_exact_.set(2, 0.0);
            }

        }

        std::string name() const override
        {
            return "Gaussaian";
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

            assert(point.size() == 3);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            result = 0.0;
            for(SizeType i = 1; i <=15; i++)
            {
                Scalar a = ( 3.5 - (0.5 * ( i - 1.0 )) - z );
                Scalar b = (x* std::exp(- 0.5 * y * a*a )) - p(i);
                result += b*b;
            }

            return true;
        }

        bool gradient(const Vector &point, Vector &g) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 3);
            g = zeros(3);

            const Read<Vector> read(point);
            const Write<Vector> write(g);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            Scalar a = 0.0;
            Scalar b = 0.0;
            Scalar c = 0.0;

            for(SizeType i =1; i <=15; i++)
            {
                Scalar d1 = 0.5 * ( i - 1.);
                Scalar d2 = 3.5 - d1 - z;
                Scalar arg = - 0.5 * y * d2 * d2;
                Scalar t = x * std::exp(arg) - p(i);

                a += 2.0 * std::exp(arg) * t;
                b -= x * std::exp(arg) * t * d2 * d2;
                c += 2.0 * x * y * std::exp(arg) * t * d2;
            }

            g.set(0, a);
            g.set(1, b);
            g.set(2, c);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 3);
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

            for(SizeType i =1; i <=15; i++)
            {
                Scalar d1 = 0.5 * ( i - 1.0);
                Scalar d2 = 3.5 - d1 - z;
                Scalar arg = 0.5 * y * d2 * d2;
                Scalar r = std::exp ( - arg );
                Scalar t = (x * r) - p(i);
                Scalar t1 = (2.0 * x * r) - p(i);

                term11 += r * r;
                term22 += r * t1 * d2*d2*d2*d2;
                term33 += r * ((y * t1 * d2 * d2) - t);
                term21 -= r * t1 * d2 * d2;
                term31 += d2 * r * t1;
                term32 += d2 * r * ( t - (arg * t1));
            }

            term11 *= 2.0;
            term22 *= 0.5 * x;
            term33 *= 2.0 * x * y;
            term31 *= 2.0 * y;
            term32 *= 2.0 * x;

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
            return 1.12793e-8;
        }


        private:
            Scalar p(const SizeType & i) const
            {
                if(i==1 || i ==15)
                    return 0.0009;
                else if(i==2 || i == 14)
                    return 0.0044;
                else if(i==3 || i == 13)
                    return 0.0175;
                else if(i==4 || i == 12)
                    return 0.0540;
                else if(i==5 || i == 11)
                    return 0.1295;
                else if(i==6 || i == 10)
                    return 0.2420;
                else if(i==7 || i == 9)
                    return 0.3521;
                else if(i==8)
                    return 0.3989;
                else
                {
                    utopia_error("Gaussian09::p():: we should never reach here... \n");
                    return 0;
                }
            }

    private:
        Vector x_init_;
        Vector x_exact_;

    };
}

#endif //UTOPIA_GAUSSIAN_09