#ifndef UTOPIA_SOLVER_PENALTY1_23
#define UTOPIA_SOLVER_PENALTY1_23

#include <vector>
#include <assert.h>
#include "utopia_Function.hpp"


namespace utopia 
{
    template<class Matrix, class Vector>
    class PenaltyI23 final: public UnconstrainedTestFunction<Matrix, Vector> 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        PenaltyI23(const SizeType & n_loc): n_loc_(n_loc) 
        {
            x_init_ = local_zeros(n_loc_);
            x_exact_ = local_values(n_loc_, 0.15812); // depends on size.. this is valid for n=10
            
            {
                const Write<Vector> write1(x_init_);
                each_write(x_init_, [](const SizeType i) -> double 
                { 
                    return i+1; 
                }   );
                 
            }
        }

        std::string name() const override
        {
            return "Penalty I"; 
        }

        SizeType dim() const override
        {
            return n_loc_; 
        }

        bool parallel() const override
        {
            return true; 
        }

        bool exact_sol_known() const override
        {
            return false; 
        }


        bool value(const Vector &x, Scalar &result) const override 
        {
            assert(local_size(x).get(0) == this->dim());
            
            Scalar alpha = 0.00001;
            Scalar t1 = -0.25 + dot(x,x); 

            Vector help = x - local_values(local_size(x).get(0), 1.0); 
            Scalar t2 = dot(help, help); 

            result = (alpha * t2) + (t1*t1); 

            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override 
        {
            assert(local_size(x).get(0) == this->dim());
            
            Scalar alpha = 0.00001;
            Scalar t1 = -0.25 + dot(x,x); 
            Vector help = x - local_values(local_size(x).get(0), 1.0); 

            g = 2.0* alpha * help; 
            g += 4.0 * t1 * x; 

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override 
        {
            assert(local_size(x).get(0) == this->dim());
            
            H = outer(x,x); 
            H *= 8.0; 

            {
                const Read<Vector> read(x);
                const Write<Matrix> write(H);

                Scalar alpha = 0.00001;
                Scalar t1 = -0.25 + dot(x,x); 
                Scalar d = (2.0 * alpha) + (4.0  *t1); 

                auto r = row_range(H); 
                for(auto i = r.begin(); i != r.end(); ++i)
                {
                    H.set(i,i, d + 8.0 * x.get(i) * x.get(i)); 
                }                
            }

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
            return 7.08765e-5; // if n=10
        }

    private: 
        SizeType n_loc_; 
        Vector x_init_; 
        Vector x_exact_; 

    };    
    
}

#endif //UTOPIA_SOLVER_PENALTY1_23
