#ifndef UTOPIA_SOLVER_BIGGS_18
#define UTOPIA_SOLVER_BIGGS_18

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia
{

    template<class Matrix, class Vector>
    class Biggs18 final: public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        using Traits   = utopia::Traits<Vector>;
        using Scalar   = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm     = typename Traits::Communicator;

        Biggs18()
        {

            x_init_.zeros(serial_layout(dim()));
            x_exact_.zeros(serial_layout(dim()));

            {
                const Write<Vector> write1(x_init_);
                const Write<Vector> write2(x_exact_);

                x_init_.set(0, 1.0);
                x_init_.set(1, 2.0);
                x_init_.set(2, 1.0);
                x_init_.set(3, 1.0);
                x_init_.set(4, 1.0);
                x_init_.set(5, 1.0);

                // x_init_.set(0, 1.0);
                // x_init_.set(1, 9.0);
                // x_init_.set(2, 1.0);
                // x_init_.set(3, 4.0);
                // x_init_.set(4, 3.0);
                // x_init_.set(5, 3.0);

                x_exact_.set(0, 1.0);
                x_exact_.set(1, 10.0);
                x_exact_.set(2, 1.0);
                x_exact_.set(3, 5.0);
                x_exact_.set(4, 4.0);
                x_exact_.set(5, 3.0);
            }

        }

        std::string name() const override
        {
            return "Biggs Exp6";
        }

        SizeType dim() const override
        {
            return 6.0;
        }

        bool exact_sol_known() const override
        {
            return false;
        }



        bool value(const Vector &point, typename Vector::Scalar &result) const override
        {
           if( point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 6);

            {
                const Read<Vector> read(point);

                const Scalar x1 = point.get(0);
                const Scalar x2 = point.get(1);
                const Scalar x3 = point.get(2);
                const Scalar x4 = point.get(3);
                const Scalar x5 = point.get(4);
                const Scalar x6 = point.get(5);


                result = 0.0;
                for(SizeType i = 1; i <=13; i++)
                {
                    Scalar c = -i*0.1;
                    Scalar b = (x3 * std::exp(c* x1)) - (x4* std::exp(c*x2))+ (x6 * std::exp(c*x5)) - std::exp (c)+ (5.0* std::exp(10.0*c)) - (3.0 * std::exp(4.0*c));
                    result += b*b;
                }

            }

            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override
        {
           if( point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 6);

            if(empty(result)){
                result = zeros(layout(point));
            }

            {
                const Read<Vector> read(point);
                const Write<Vector> write(result);

                const Scalar x1 = point.get(0);
                const Scalar x2 = point.get(1);
                const Scalar x3 = point.get(2);
                const Scalar x4 = point.get(3);
                const Scalar x5 = point.get(4);
                const Scalar x6 = point.get(5);

                Scalar g1=0, g2=0, g3=0, g4=0, g5=0, g6 =0;

                for(SizeType i = 1; i <=13; i++)
                {
                    Scalar c = -i*0.1;
                    Scalar fi = (x3 * std::exp(c* x1)) - (x4* std::exp(c*x2))+ (x6 * std::exp(c*x5)) - std::exp (c)+ (5.0* std::exp(10.0*c)) - (3.0 * std::exp(4.0*c));

                    Scalar dfdx1 =   c*x3* std::exp(c*x1);
                    Scalar dfdx2 = - c*x4* std::exp(c*x2);
                    Scalar dfdx3 =         std::exp(c*x1);
                    Scalar dfdx4 = -       std::exp(c*x2);
                    Scalar dfdx5 =   c*x6* std::exp(c*x5);
                    Scalar dfdx6 =         std::exp(c*x5);

                    g1+= 2.0 * fi * dfdx1;
                    g2+= 2.0 * fi * dfdx2;
                    g3+= 2.0 * fi * dfdx3;
                    g4+= 2.0 * fi * dfdx4;
                    g5+= 2.0 * fi * dfdx5;
                    g6+= 2.0 * fi * dfdx6;
                }

                result.set(0, g1);
                result.set(1, g2);
                result.set(2, g3);
                result.set(3, g4);
                result.set(4, g5);
                result.set(5, g6);
            }
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override
        {
           if( point.comm().size() > 1){
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            assert(point.size() == 6);
            // result.dense(serial_layout(6, 6));
            if(empty(result)){
                result.dense(serial_layout(6, 6));
            }

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x1 = point.get(0);
            const Scalar x2 = point.get(1);
            const Scalar x3 = point.get(2);
            const Scalar x4 = point.get(3);
            const Scalar x5 = point.get(4);
            const Scalar x6 = point.get(5);

            Scalar term11= 0, term12= 0, term13= 0, term14= 0, term15= 0, term16 = 0;
            Scalar term22= 0, term23= 0, term24= 0, term25= 0, term26 = 0;
            Scalar term33= 0, term34= 0, term35= 0, term36 = 0;
            Scalar term44= 0, term45= 0, term46 = 0, term55= 0, term56 = 0, term66 = 0;


            for(SizeType i = 1; i <=13; i++)
            {
                Scalar c = -i*0.1;
                Scalar fi = (x3 * std::exp(c* x1)) - (x4* std::exp(c*x2))+ (x6 * std::exp(c*x5)) - std::exp (c)+ (5.0* std::exp(10.0*c)) - (3.0 * std::exp(4.0*c));

                Scalar dfdx1 = c*x3*std::exp(c*x1);
                Scalar d2fdx11 = c*c*x3* std::exp(c*x1);
                Scalar d2fdx13 = c*std::exp(c*x1);
                Scalar dfdx2 = -c*x4*std::exp(c*x2);
                Scalar d2fdx22 = -c*c*x4* std::exp(c*x2);
                Scalar d2fdx24 = -c*std::exp(c*x2);
                Scalar dfdx3 = std::exp(c*x1);
                Scalar dfdx4 = -std::exp(c*x2);
                Scalar dfdx5 = c*x6*std::exp(c*x5);
                Scalar d2fdx55 = c*c*x6* std::exp(c*x5);
                Scalar d2fdx56 =c*std::exp(c*x5);
                Scalar dfdx6 = std::exp(c*x5);

                term11 += (2.0 * dfdx1 * dfdx1) + (2.0 * fi * d2fdx11);
                term12 += 2.0 * dfdx2 * dfdx1;
                term13 += (2.0 * dfdx3 * dfdx1) + (2.0 * fi * d2fdx13);
                term14 += 2.0 * dfdx4 * dfdx1;
                term15 += 2.0 * dfdx5 * dfdx1;
                term16 += 2.0 * dfdx6 * dfdx1;

                term22 += (2.0 * dfdx2 * dfdx2) + (2.0 * fi * d2fdx22);
                term23 += 2.0 * dfdx3 * dfdx2;
                term24 += (2.0 * dfdx4 * dfdx2) + (2.0 * fi * d2fdx24);
                term25 += 2.0 * dfdx5 * dfdx2;
                term26 += 2.0 * dfdx6 * dfdx2;

                term33 += 2.0 * dfdx3 * dfdx3;
                term34 += 2.0 * dfdx4 * dfdx3;
                term35 += 2.0 * dfdx5 * dfdx3;
                term36 += 2.0 * dfdx6 * dfdx3;

                term44 += 2.0 * dfdx4 * dfdx4;
                term45 += 2.0 * dfdx5 * dfdx4;
                term46 += 2.0 * dfdx6 * dfdx4;

                term55 += (2.0 * dfdx5 * dfdx5) + (2.0 * fi * d2fdx55);
                term56 += (2.0 * dfdx6 * dfdx5) + (2.0 * fi * d2fdx56);

                term66 += 2.0 * dfdx6 * dfdx6;
            }

            result.set(0, 0, term11);
            result.set(0, 1, term12);
            result.set(0, 2, term13);
            result.set(0, 3, term14);
            result.set(0, 4, term15);
            result.set(0, 5, term16);

            result.set(1, 0, term12);
            result.set(1, 1, term22);
            result.set(1, 2, term23);
            result.set(1, 3, term24);
            result.set(1, 4, term25);
            result.set(1, 5, term26);

            result.set(2, 0, term13);
            result.set(2, 1, term23);
            result.set(2, 2, term33);
            result.set(2, 3, term34);
            result.set(2, 4, term35);
            result.set(2, 5, term36);

            result.set(3, 0, term14);
            result.set(3, 1, term24);
            result.set(3, 2, term34);
            result.set(3, 3, term44);
            result.set(3, 4, term45);
            result.set(3, 5, term46);

            result.set(4, 0, term15);
            result.set(4, 1, term25);
            result.set(4, 2, term35);
            result.set(4, 3, term45);
            result.set(4, 4, term55);
            result.set(4, 5, term56);

            result.set(5, 0, term16);
            result.set(5, 1, term26);
            result.set(5, 2, term36);
            result.set(5, 3, term46);
            result.set(5, 4, term56);
            result.set(5, 5, term66);


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
            return 5.65565;
        }


    private:
        Vector x_init_;
        Vector x_exact_;


    };



}

#endif //UTOPIA_SOLVER_BIGGS_18
