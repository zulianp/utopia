#ifndef UTOPIA_BROWN_DENNIS_16
#define UTOPIA_BROWN_DENNIS_16

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class BrownDennis16 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        BrownDennis16()
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(4);
            x_exact_ = zeros(4);

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                // x_init_.set(0, -11.0);
                // x_init_.set(1, 13.0);
                // x_init_.set(2, -0.5);
                // x_init_.set(3, 0.2);

                x_init_.set(0, 25.0);
                x_init_.set(1, 5.0);
                x_init_.set(2, -5.0);
                x_init_.set(3, -1.0);

                x_exact_.set(0, -11.5844);
                x_exact_.set(1, 13.1999);
                x_exact_.set(2, -0.406200);
                x_exact_.set(3, 0.240998);
            }

        }

        std::string name() const override
        {
            return "Brown and Dennis";
        }

        SizeType dim() const override
        {
            return 4;
        }


        bool exact_sol_known() const override
        {
            return false;  // just because we can not fit into precision
        }


        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 4);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);
            const Scalar w = point.get(3);

            result = 0.0;
            for(SizeType i = 1; i <=20; i++)
            {
                Scalar c = i*0.2;
                Scalar f1 = x + (c * y) - std::exp(c);
                Scalar f2 = z + (std::sin(c) * w) - std::cos(c);

                result += std::pow(f1, 4.) + (2.0 * f1*f1 * f2*f2) + std::pow(f2, 4.);
            }

            return true;
        }

        bool gradient(const Vector &point, Vector &g) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 4);
            g = zeros(4);

            const Read<Vector> read(point);
            const Write<Vector> write(g);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);
            const Scalar w = point.get(3);

            Scalar g1 = 0.0;
            Scalar g2 = 0.0;
            Scalar g3 = 0.0;
            Scalar g4 = 0.0;

            for(SizeType i =1; i <=20; i++)
            {
                Scalar c = i*0.2;
                Scalar f1 = x + (c * y) - std::exp(c);
                Scalar f2 = z + (std::sin(c) * w) - std::cos(c);

                Scalar  df2dx4 = std::sin (c);

                g1 += 4.0 * ((std::pow(f1,3)) + (f1 * f2 * f2));
                g2 += 4.0 * ((std::pow(f1,3)* c) + (f1 * f2 * f2 * c));
                g3 += 4.0 * ((f1 * f1 * f2 ) + (std::pow(f2,3)));
                g4 += 4.0 * ((f1 * f1 * f2 * df2dx4) + (std::pow(f2,3) * df2dx4));
            }

            g.set(0, g1);
            g.set(1, g2);
            g.set(2, g3);
            g.set(3, g4);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
            if( mpi_world_size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size().get(0) == 4);
            result = zeros(4,4);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);
            const Scalar w = point.get(3);

            Scalar term11 = 0.0;
            Scalar term22 = 0.0;
            Scalar term33 = 0.0;
            Scalar term21 = 0.0;
            Scalar term31 = 0.0;
            Scalar term32 = 0.0;
            Scalar term41 = 0.0;
            Scalar term42 = 0.0;
            Scalar term34 = 0.0;
            Scalar term44 = 0.0;

            for(SizeType i =1; i <=20; i++)
            {
                Scalar c = i*0.2;
                Scalar f1 = x + (c * y) - std::exp(c);
                Scalar f2 = z + (std::sin(c) * w) - std::cos(c);

                Scalar  df2dx4 = std::sin (c);

                term11 += (12.0 * f1 * f1) +  (4.0 * f2 * f2);
                term22 += (12.0 * f1 * f1 * c * c) +  (4.0 * f2 * f2 * c);
                term33 += (4.0 * f1 * f1) +  (12.0 * f2 * f2);
                term21 += (12.0 * f1 * f1 * c) +  (4.0 * f2 * f2 * c);
                term31 += 8.0 * f1 * f2;
                term32 += 8.0 * f1 * f2 * c;
                term41 += 8.0 * f1 * f2 * df2dx4;
                term42 += 8.0 * f1 * f2 * c * df2dx4;
                term34 += (4.0 * f1 * f1 * df2dx4) + (12.0 * f2 * f2 * df2dx4);
                term44 += (4.0 * f1 * f1 * df2dx4 * df2dx4) + (12.0 * f2 * f2 * df2dx4 * df2dx4);
            }

            result.set(0, 0, term11);
            result.set(0, 1, term21);
            result.set(0, 2, term31);
            result.set(0, 3, term41);

            result.set(1, 0, term21);
            result.set(1, 1, term22);
            result.set(1, 2, term32);
            result.set(1, 3, term42);

            result.set(2, 0, term31);
            result.set(2, 1, term32);
            result.set(2, 2, term33);
            result.set(2, 3, term34);

            result.set(3, 0, term41);
            result.set(3, 1, term42);
            result.set(3, 2, term34);
            result.set(3, 3, term44);

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
            return 85822.2;
        }

    private:
        Vector x_init_;
        Vector x_exact_;

    };
}

#endif //UTOPIA_BROWN_DENNIS_16