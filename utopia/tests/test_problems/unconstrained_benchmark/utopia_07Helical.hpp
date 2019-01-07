#ifndef UTOPIA_HELICAL_07
#define UTOPIA_HELICAL_07

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_UnconstrainedTestFunction.hpp"


namespace utopia
{   
    template<class Matrix, class Vector>
    class Hellical07 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        Hellical07() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(3);
            x_exact_ = zeros(3);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                x_init_.set(0, -1.0);
                x_init_.set(1, 0.0);
                x_init_.set(2, 0.0);

                x_exact_.set(0, 1.0);
                x_exact_.set(1, 0.0);                
                x_exact_.set(2, 0.0);                                
            }

        }

        std::string name() const override
        {
            return "Hellical valley"; 
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

            Scalar a = 10.0 * (z - 10.0 * theta(x,y));
            Scalar b = 10.0 * ( std::sqrt(x*x + y*y) - 1.0); 

            result = a*a + b*b + z*z; 
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            assert(point.size().get(0) == 3);
            result = zeros(3);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar z = point.get(2);

            const Scalar xx = x * x;
            const Scalar yy = y * y; 
            const Scalar r = std::sqrt(xx + yy); 
            const Scalar t =  z - 10.0 * theta(x,y); 
            const Scalar s1 = 5.0 * t / ( pi() * r * r );

            const Scalar a = 200.0 * ( x - (x / r) + y * s1 );
            const Scalar b = 200.0 * ( y - (y / r) - x * s1 );
            const Scalar c = 2.0 * ( 100.0 * t + z);

            result.set(0, a);
            result.set(1, b);
            result.set(2, c);

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

            const Scalar xx = x * x;
            const Scalar yy = y * y; 
            const Scalar xy = x * y; 
            const Scalar xxyy = xx + yy; 

            const Scalar pixy = pi() * ( xx + yy ); 
            const Scalar xxyy32 = std::pow(xxyy , 3./2.); 


            const Scalar th = theta(x,y); 
            Scalar h1 = pi() * (xxyy);     
            Scalar h2 = h1 * (xxyy);     

            Scalar  term11 = 200.0 - 200.0 * yy * ( 1.0 / xxyy32 - 25.0 / ( h1 *h1 ));
                    term11 -= 2000.0 * xy * ( z - 10.0 * th )/h2;

            const Scalar mixed23 = - 1000.0 * x / pixy;
            const Scalar mixed13 = 1000.0 * y / pixy;


            Scalar  term12 = 200.0 * xy / xxyy32;
                    term12 += 1000.0 /h2; 
                    term12 *= ( ( z - 10.0 * th ) * ( xx - yy ) - 5.0 * xy / pi() );


            Scalar  term22 = 200.0 - 200.0 * xx * ( 1.0 / xxyy32 - 25.0 / (h1*h1));
                    term22 += 2000.0 * xy * ( z - 10.0 * th )/h2;


            result.set(0, 0, term11);
            result.set(0, 1, term12);
            result.set(0, 2, mixed13);

            result.set(1, 0, term12);
            result.set(1, 1, term22);
            result.set(1, 2, mixed23);
            
            result.set(2, 0, mixed13);
            result.set(2, 1, mixed23);
            result.set(2, 2, 202.0);    

            
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
            Scalar theta(const Scalar & x1, const Scalar & x2) const
            {
              if ( 0.0 < x1 )
                return 0.5 * std::atan ( x2 / x1 ) / pi();
              else if ( x1 < 0.0 )
                return 0.5 * std::atan ( x2 / x1 ) / pi() + 0.5;
              else if ( 0.0 < x2 )
                return 0.25;
              else if ( x2 < 0.0 )
                return - 0.25;
              else
                return 0.0;
            }

            constexpr Scalar pi() const
            { 
                return 3.141592653589793238462643383279502884; 
            }


    private: 
        Vector x_init_; 
        Vector x_exact_; 

    };
}

#endif //UTOPIA_HELICAL_07