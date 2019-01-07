#ifndef UTOPIA_BOX_12
#define UTOPIA_BOX_12

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_UnconstrainedTestFunction.hpp"


namespace utopia
{   
    template<class Matrix, class Vector>
    class Box12 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Box12() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(3);
            x_exact_ = zeros(3);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                // trust region does not converge to global sol, starting here...
                x_init_.set(0, 1.0);
                x_init_.set(1, 10.0);
                x_init_.set(2, 5.0);

                // x_init_.set(0, 1.0);
                // x_init_.set(1, 9.5);
                // x_init_.set(2, 1.0);                

                x_exact_.set(0, 1.0);
                x_exact_.set(1, 10.0);                
                x_exact_.set(2, 1.0);           
                // also {10, 1, -1} or {x1, x1, 0}                     
            }

        }

        std::string name() const override
        {
            return "Box three-dimensional"; 
        }

        SizeType dim() const override
        {
            return 3; 
        }


        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 3);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            result = 0.0; 
            for(SizeType i = 1; i <=10; i++)
            {   
                Scalar c = -i*0.1; 
                Scalar b = std::exp(c*x) - std::exp(c*y) - (z*( std::exp(c)- std::exp(10.0 * c)));
                result += b*b;
            }
    
            return true;
        }

        bool gradient(const Vector &point, Vector &g) const override 
        {
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

            for(SizeType i =1; i <=10; i++)
            {   
                Scalar c = -i*0.1; 
                Scalar fi = std::exp(c*x) - std::exp(c*y) - (z * (std::exp(c)- std::exp(10.0*c)));

                Scalar  dfidx1 = c*std::exp(c*x);
                Scalar  dfidx2 = -c*std::exp(c*y);
                Scalar  dfidx3 = -(std::exp(c)-std::exp(10.0*c));

                g1 += 2.0 * fi * dfidx1;
                g2 += 2.0 * fi * dfidx2;
                g3 += 2.0 * fi * dfidx3;
            }

            g.set(0, g1);
            g.set(1, g2);
            g.set(2, g3);

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
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

            for(SizeType i =1; i <=10; i++)
            {   
                Scalar c = -i*0.1; 
                Scalar fi = std::exp(c*x) - std::exp(c*y) - (z * (std::exp(c)- std::exp(10.0*c)));

                Scalar  dfidx1      = c * std::exp(c*x);
                Scalar  d2fidx11    = c * c * std::exp(c*x);
                Scalar  dfidx2      = - c * std::exp(c*y);
                Scalar  d2fidx22    = - c * c * std::exp(c*y);
                Scalar  dfidx3      = - ( std::exp(c) - std::exp(10.0*c));

                term11 += (2.0 * dfidx1 * dfidx1) + (2.0 * fi * d2fidx11);
                term22 += (2.0 * dfidx2 * dfidx2) + (2.0 * fi * d2fidx22);
                term33 += 2.0 * dfidx3 * dfidx3;
                term21 += 2.0 * dfidx1 * dfidx2;
                term31 += 2.0 * dfidx1 * dfidx3;
                term32 += 2.0 * dfidx3 * dfidx2;
            }

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

#endif //UTOPIA_BOX_12