#ifndef UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
#define UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"


namespace utopia
{

    // This function is used for implementation of test functions for unconstrained nonlinear benchmark
    // see: More, Garbow, Hillstrom - Testing unconstrained optimization software
    template<class Matrix, class Vector>
    class UnconstrainedTestFunction : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        virtual ~UnconstrainedTestFunction() { }

        virtual Vector initial_guess() const = 0;
        virtual const Vector & exact_sol() const = 0;
        virtual Scalar min_function_value() const = 0;
    };


    /**
     * @brief      Rosenbrock 2D banana function. \n 
     *             The floor of the valley follows approximately the parabola \f$ y = x^2 + 1/200 \f$.   
     *             The covariance matrix is not positive-definite. On the dashed line it is singular. 
     *             Stepping method tend to perform at least as well as gradient methods for this function.
     *
     */
    template<class Matrix, class Vector>
    class Rosenbrock final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        Rosenbrock() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(2);
            x_exact_ = zeros(2);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                x_init_.set(0, -1.2);
                x_init_.set(1, 1.0);

                x_exact_.set(0, 1.0);
                x_exact_.set(1, 1.0);                
            }

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result = 1 + 100.0 * pow(x * x - y , 2.0) + pow(x - 1 , 2.0);
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

            result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
            result.set(1, 200.0 * (y - x * x));
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
            return 0; 
        }


    private: 
        Vector x_init_; 
        Vector x_exact_; 

    };



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


    template<class Matrix, class Vector>
    class Brown04 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        Brown04() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(2);
            x_exact_ = zeros(2);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                x_init_.set(0, 1.0);
                x_init_.set(1, 1.0);

                x_exact_.set(0, 1e6);
                x_exact_.set(1, 2e-6);                
            }

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            Scalar a = x - 1e6; 
            Scalar b = y - 2e-6; 
            Scalar c = x*y - 2.0; 

            result = a*a + b*b + c*c; 
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

            Scalar a = (2.0*x* (y*y + 1)) - (4.0 * (y + 500000)); 
            Scalar b = (2.0*x* (x*y -2.0)) + 2.0*y - 2.0/500000; 

            result.set(0, a);
            result.set(1, b);

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

            const Scalar mixed = (4.0 * x * y) - 4.0; 
            const Scalar a = 2.0* (y*y +1); 
            const Scalar b = 2.0* (x*x +1); 

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

    template<class Matrix, class Vector>
    class Beale05 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        Beale05() 
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_init_ = zeros(2);
            x_exact_ = zeros(2);

            const Write<Vector> write1(x_init_);
            const Write<Vector> write2(x_exact_);
            {
                x_init_.set(0, 1.0);
                x_init_.set(1, 1.0);

                x_exact_.set(0, 3.0);
                x_exact_.set(1, 0.5);                
            }

        }

        bool value(const Vector &point, typename Vector::Scalar &result) const override 
        {
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


    template<class Matrix, class Vector>
    class Hellical07 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

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
#endif //UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
