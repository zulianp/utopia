/*
* @Author: alenakopanicakova
* @Date:   2016-05-10
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-03
*/

#ifndef UTOPIA_BACKTRACKING_HPP
#define UTOPIA_BACKTRACKING_HPP
#include "utopia_LS_Strategy.hpp"
#include "utopia_PrintInfo.hpp"


namespace utopia 
{

    /**
     * @brief      This class implements very simple backtracking, without any interpolation with Wolfe SDC.  
     *             It is mainly made for learning and testing purposes. \n
     *             For more serious implementation check Backtracing class. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class SimpleBacktracking : public LSStrategy<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

    public:


        SimpleBacktracking(const Parameters params = Parameters())

                        : LSStrategy<Matrix, Vector>(params)

        {
            set_parameters(params); 
        }

        /**
         * @brief      Gets the step-size alpha
         *
         * @param      fun      The nonlinear function. 
         * @param[in]  g        The gradient.
         * @param[in]  x        The current iterate.
         * @param[in]  p_k      The new step.
         * @param      alpha_k  The step-size. 
         *
         * @return     
         */

        bool get_alpha(LeastSquaresFunction<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override 
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }


        bool get_alpha(Function<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }

        template<class FunctionT>
        bool get_alpha_aux(FunctionT &fun, const Vector &g, const Vector& x, const Vector &p_k, Scalar &alpha_k)
        {
            Vector x_0 = x, x_k = x;
            Scalar E_k, E_k1, g_p;

            fun.value(x_k, E_k);
            alpha_k = 1;
            g_p =  dot(g, p_k);

            x_k = x_k + p_k;
            fun.value(x_k, E_k1);

            SizeType it = 0; 

            if(verbose_)
                PrintInfo::print_init("SIMPLE_BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| E_k1 ||"}); 

            // Wolfe conditions                        
            while( E_k1 < (E_k + c1_ * alpha_k * g_p) && it < max_it_  && alpha_k > 1e-6)
            {

                alpha_k = alpha_k * rho_;
                x_k = x_0 + alpha_k * p_k;
                fun.value(x_k, E_k1);
                it++; 
                
                if(verbose_)
                    PrintInfo::print_iter_status(it, {E_k1}); 

            }

            std::cout<<"alpha_k:  "<< alpha_k << "  \n"; 
            return true;
        }


        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        bool set_parameters(const Parameters params) override
        {            
            verbose_    = params.line_search_inner_verbose(); 
            c1_         = params.c1(); 
            max_it_     = params.n_line_search_iters(); 
            rho_        = params.ls_rho();
            alpha_min_   = params.alpha_min(); 
            
            return true; 
        } 


    private:
        bool verbose_;      /*!< Verbose inside of LS strategy.  */  
        Scalar c1_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */  
        Scalar max_it_;     /*!< Maximum of the iterations inside of LS strategy.  */  
        Scalar rho_;        /*!< Contraction factor.   */  
        Scalar alpha_min_;  /*!< Minimum allowed step-size.   */  

    };

}


#endif //UTOPIA_BACKTRACKING_HPP

