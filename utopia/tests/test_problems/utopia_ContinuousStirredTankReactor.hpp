#ifndef UTOPIA_CSR_EXAMPLE_HPP
#define UTOPIA_CSR_EXAMPLE_HPP

#include <vector>
#include <assert.h>
#include <cmath>
#include "utopia_Function.hpp"



namespace utopia
{

    /**
     * @brief       Example taken from "Numerical methods for chemical engineers with MATLAB applications" by Constantinides, Mostoufi, et al.
     *              Example used for testing ASTRUM solver in "A new adaptive approach to pseudo-transient continuation method" by Kopanicakova, Krause, Deuflhard
     *              
     *             Note: - example can run only in serial  
     *              
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class ContinuousStirredReactor : public UnconstrainedTestFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        ContinuousStirredReactor(): 
        v_(1.0),
        V_(100.0),
        k1_(1.0),
        k2_(1.0),
        u1_in_(1.0),
        u2_in_(2.0),
        u3_in_(0.0),
        u4_in_(0.0)
        {
            assert(!utopia::is_parallel<Matrix>::value || mpi_world_size() == 1 && "does not work for parallel matrices");

            x_exact_ = zeros(4);

            {
                const Write<Vector> write2(x_exact_);

                // valid only for constants choosen in the constructor
                x_exact_.set(0, 0.056613650336565);
                x_exact_.set(1, 0.166635845605265);
                x_exact_.set(2, 0.053408544932134);
                x_exact_.set(3, 0.889977804731300);
            }

        }

        std::string name() const override
        {
            return "Continuous Stirred Reactor";
        }

        SizeType dim() const override
        {
            return 4;
        }

        bool exact_sol_known() const override
        {
            return true;
        }

        Vector initial_guess() const override
        {
            return values(4, 0); 
        }

        const Vector & exact_sol() const override
        {
            return x_exact_;
        }

        Scalar min_function_value() const override
        {
            return 0.0;
        }


        bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            // merit function
            assert(x.size() == 4);
            Vector g = values(2, 0.0);
            gradient(x, g);
            result = 0.5 * norm2(g);
            return true;
        }

        bool gradient(const Vector &x, Vector &g) const override
        {
            assert(x.size() == 4);
            g = zeros(4);

            {
                const Read<Vector> read(x);
                const Write<Vector> write(g);

                const Scalar x1 = x.get(0);
                const Scalar x2 = x.get(1);
                const Scalar x3 = x.get(2);
                const Scalar x4 = x.get(3);

                g.set(0, v_* (x1 - u1_in_) + V_ * (k1_*x1*x2));
                g.set(1, v_* (x2 -u2_in_)  + V_ * ((k1_*x1*x2)+(k2_*x2*x3))); 
                g.set(2, v_ *(x3 - u3_in_) + V_ * ((-k1_*x1*x2)+(k2_*x2*x3)));
                g.set(3, v_*(x4 -u4_in_)   + V_ * (-k2_*x2*x3));
            }

            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            assert(x.size() == 4);

            H = zeros(4,4);

            {
                const Read<Vector> read(x);
                const Write<Matrix> write(H);

                const Scalar x1 = x.get(0);
                const Scalar x2 = x.get(1);
                const Scalar x3 = x.get(2);
                // const Scalar x4 = x.get(3);

                H.set(0, 0, (v_)+ (V_*k1_*x2));
                H.set(0, 1, V_*k1_*x1);
                H.set(0, 2, 0.0);
                H.set(0, 3, 0.0);

                H.set(1, 0, V_* k1_ *x2);
                H.set(1, 1, (v_)+ (V_*k1_*x1)+(V_ * k2_*x3));
                H.set(1, 2, V_*k2_*x2);
                H.set(1, 3, 0.0);

                H.set(2, 0, V_*-k1_*x2);
                H.set(2, 1, V_*((-k1_*x1)+(k2_*x3)));
                H.set(2, 2, (v_)+(V_*k2_*x2));
                H.set(2, 3, 0.0);         

                H.set(3, 0, 0.0);
                H.set(3, 1, V_*-k2_*x3);
                H.set(3, 2, V_*-k2_*x2);
                H.set(3, 3, v_);
            }

            return true;
        }


        void get_initial_guess(Vector & x, const Scalar & value = 0)
        {
            x = values(4, value);
        }


        void vol_flow_rate(const Scalar & v)
        {
            v_ = v;
        }

        void reactor_volume(const Scalar & V)
        {
            V_ = V;
        }     

        void rate_const_1(const Scalar & k1)
        {
            k1_ = k1;
        }     

        void rate_const_2(const Scalar & k2)
        {
            k2_ = k2;
        } 

        void inlet_concentration_A(const Scalar & u1_in)
        {
            u1_in_ = u1_in;
        }    

        void inlet_concentration_B(const Scalar & u2_in)
        {
            u2_in_ = u2_in;
        }                 

        void inlet_concentration_C(const Scalar & u3_in)
        {
            u3_in_ = u3_in;
        }    

        void inlet_concentration_D(const Scalar & u4_in)
        {
            u4_in_ = u4_in;
        }            


        Scalar vol_flow_rate() const
        {
            return v_;
        }

        Scalar reactor_volume() const
        {
            return V_;
        }     

        Scalar rate_const_1() const
        {
            return  k1_;
        }     

        Scalar rate_const_2() const
        {
            return  k2_;
        } 

        Scalar inlet_concentration_A() const 
        {
            return  u1_in_;
        }    

        Scalar inlet_concentration_B() const
        {
            return u2_in_;
        }                 

        Scalar inlet_concentration_C() const 
        {
            return  u3_in_;
        }    

        Scalar inlet_concentration_D() const 
        {
            return u4_in_;
        }  



    private: 
        Scalar v_; 
        Scalar V_; 
        Scalar k1_; 
        Scalar k2_; 
        Scalar u1_in_; 
        Scalar u2_in_; 
        Scalar u3_in_; 
        Scalar u4_in_; 

        Vector x_exact_; 

    };

}

#endif //UTOPIA_CSR_EXAMPLE_HPP
