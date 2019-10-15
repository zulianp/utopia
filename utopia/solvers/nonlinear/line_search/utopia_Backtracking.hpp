#ifndef UTOPIA_QUAD_CUB_BACKTRACKING_HPP
#define UTOPIA_QUAD_CUB_BACKTRACKING_HPP

#include "utopia_LS_Strategy.hpp"
#include "utopia_Function.hpp"
#include "utopia_FunctionNormalEq.hpp"

namespace utopia {

    /**
     * @brief

        This is the algorithm described in the book
        Dennis and Schnabel "Numerical Methods for Nonlinear Equations
        and Unconstrained Optimization", Prentice-Hall (1983),
        reprinte by SIAM (1996), Section 6.3.2.
        @todo check params naming properly...
     */
    template<class Vector, int Backend = Traits<Vector>::Backend>
    class Backtracking final : public LSStrategy<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;


    public:
        Backtracking(): LSStrategy<Vector>(), c2_(1e-8)

        {

        }

        void read(Input &in) override
        {
            LSStrategy<Vector>::read(in);
            in.get("c2", c2_);
        }

        void print_usage(std::ostream &os) const override
        {
            LSStrategy<Vector>::print_usage(os);
            this->print_param_usage(os, "c2", "double", "Constant used for Wolfe conditions.", "1e-8");
        }

        void c2(const Scalar  & c)
        {
            c2_ = c;
        }


        Scalar c2() const
        {
            return c2_;
        }


        /**
         * @brief      Get the alpha_k on given iterate. We are using quadratic and qubic interpolation as part of backtracking.
         *             For checking decrease conditions, we are using Wolfe conditions.
         *
         * @param      fun      The fun with eval.
         * @param[in]  g        The gradient.
         * @param[in]  x        The current iterate.
         * @param[in]  d        The descent direction/ usually Newton step.
         * @param      alpha_k  The new step size
         *
         * @return
         */

        bool get_alpha(LeastSquaresFunctionBase<Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }


        bool get_alpha(FunctionBase<Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }

        template<class FunctionT>
        bool get_alpha_aux_home_made(FunctionT &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha)
        {
            Scalar alpha_c, alpha_p, dg = dot(d,g);
            Scalar f, f0, fc, fp, t1, t2, t3, a, b, disc;
            alpha = 1.0;

            if(dg >= 0.0)
            {
                if(mpi_world_rank() == 0) {
                    std::cerr<< "utopia::LS::backtracking:: d is not descent direction \n";
                    assert(false);
                }

                alpha = 0.0;
                return false;
            }

            fun.value(x, f);
            f0 = f;
            fc = f;
            alpha_c = alpha;

            Scalar it = 0;

            if(this->verbose())
                PrintInfo::print_init("BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| alpha ||"});

            while(alpha > c2_ && it < this->max_it())
            {
                x_k = x + alpha * d;
                fun.value(x_k, f);

                // check decrease condition (wolfe condition)
                if(f < f0 + this->c1() * alpha * dg )
                {

                    return true;
                }

                alpha_p = alpha_c;
                alpha_c = alpha;
                fp = fc;
                fc = f;

                //  compute next step size alpha
                if(it == 0)
                {
                    alpha = - dg / (2 * (fc - f0 -dg));
                    it++;
                }
                else
                {
                    // all subsequent backtracks: cubic fit
                    t1 = fc - f0 - alpha_c * dg;
                    t2 = fp - f0 - alpha_p * dg;
                    t3 = 1 / (alpha_c - alpha_p);

                    a = t3 * ( t1/std::pow(alpha_c, 2) - t2/std::pow(alpha_p, 2) );
                    b = t3 * ( t2 * alpha_c/ std::pow(alpha_p, 2) - t1 * alpha_p/std::pow(alpha_c, 2) );
                    disc = std::pow(b, 2)  - 3 * a * dg;

                    if( a != 0)
                    {
                        // cubic has unique minimum
                        alpha = (-b + std::sqrt( disc ) ) / (3 * a );
                    }
                    else
                    {
                        // cubic is a quadratic
                        alpha = -dg / ( 2 * b );
                    }
                }

                //  saveguard the step size
                if(alpha > 0.5 * alpha_c)
                {
                    alpha = 0.5 * alpha_c;
                }

                if(alpha < alpha_c/10)
                {
                    alpha = alpha_c/10;
                }
                it++;
                if(this->verbose())
                    PrintInfo::print_iter_status({it, alpha});
            }

            return true;
        }


    private:
        Scalar c2_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */
        Vector x_k; 

    };
}


#endif //UTOPIA_QUAD_CUB_BACKTRACKING_HPP

