/*
* @Author: alenakopanicakova
* @Date:   2016-06-01
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-22
*/

#ifndef UTOPIA_QUAD_CUB_BACKTRACKING_HPP
#define UTOPIA_QUAD_CUB_BACKTRACKING_HPP
#include "utopia_LS_Strategy.hpp"


#ifdef WITH_PETSC
    #include <petscsnes.h>
    #include "utopia_PETScFunction.hpp"
#endif //WITH_PETSC 

namespace utopia 
{

    /**
     * @brief      

        This is the algorithm described in the book
        Dennis and Schnabel "Numerical Methods for Nonlinear Equations 
        and Unconstrained Optimization", Prentice-Hall (1983), 
        reprinte by SIAM (1996), Section 6.3.2.
        @todo check params naming properly... 
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Backtracking : public LSStrategy<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;


    public:


        Backtracking(const Parameters params = Parameters() )
                        : LSStrategy<Matrix, Vector>(params)

        {
            set_parameters(params); 
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

        bool get_alpha(LeastSquaresFunction<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override 
        {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }


        bool get_alpha(Function<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux_home_made(fun, g, x, d, alpha);
        }

        template<class FunctionT>
        bool get_alpha_aux_home_made(FunctionT &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) 
        {
            Vector x_0 = x, x_k = x;
            Scalar alpha_c, alpha_p, dg = dot(d,g);
            Scalar f, f0, fc, fp, t1, t2, t3, a, b, disc; 
            alpha = 1; 

            if(dg >= 0)
            {
                if(mpi_rank == 0) {
                    std::cerr<< "utopia::LS::backtracking:: d is not descent direction \n"; 
                    assert(false);
                }

                alpha = 0.0;
                return false;
            }

            fun.value(x_k, f);
            f0 = f;
            fc = f;
            alpha_c = alpha; 

            Scalar it = 0; 

            if(verbose_)
                PrintInfo::print_init("BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| alpha ||"}); 

            while(alpha > c2_ && it < max_it_)
            {
                x_k = x_0 + alpha * d;
                fun.value(x_k, f);

                // check decrease condition (wolfe condition)
                if(f < f0 + c1_ * alpha * dg )
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
                if(verbose_)
                    PrintInfo::print_iter_status({it, alpha}); 
            }

            return true;
        }


    bool set_parameters(const Parameters params) override
    {
        verbose_    = params.line_search_inner_verbose(); 
        c1_         = params.c1(); 
        c2_         = params.c2(); 
        max_it_     = params.n_line_search_iters(); 
        rho_        = params.ls_rho();
        alpha_min_   = params.alpha_min(); 
        
        return true; 
    } 

    private:
        SizeType mpi_size = mpi_world_size();
        SizeType mpi_rank = mpi_world_rank();

        bool verbose_;      /*!< Verbose inside of LS strategy.  */  
        Scalar c1_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */  
        Scalar c2_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */  
        Scalar max_it_;     /*!< Maximum of the iterations inside of LS strategy.  */  
        Scalar rho_;        /*!< Contraction factor.   */  
        Scalar alpha_min_;  /*!< Minimum allowed step-size.   */  

    };


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////    PETSC VERSION  ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_PETSC
    template<typename Matrix, typename Vector>
    class Backtracking<Matrix, Vector, PETSC> : public LSStrategy<Matrix, Vector>
    {

        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;

        
        bool get_alpha(LeastSquaresFunction<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override 
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }


        bool get_alpha(Function<Matrix, Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }


        template<class FunctionT>
        bool get_alpha_aux(FunctionT &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) 
        {

            SNESLineSearch      linesearch;
            if(dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun))
            {
                PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);

                SNES snes; 
                fun_petsc->getSNES(snes); 

                SNESSetType(snes, SNESNEWTONLS); 
                SNESGetLineSearch(snes, &linesearch);
                SNESLineSearchSetType(linesearch,  SNESLINESEARCHBT); 
                SNESLineSearchSetOrder(linesearch, SNES_LINESEARCH_ORDER_CUBIC); 

                PetscReal fnorm; 
                PetscReal lambda; 
                VecNorm(raw_type(g), NORM_2, &fnorm); 


                // if we reset objective, computation is done in Least-squares manner 
                // SNESSetObjective(snes, NULL, NULL);  

                // we do this, since current solution and F are going to be changed in apply routine...
                // is it more efficient, but current solvers are set up differently
                Vector x_p = x; 
                Vector g_p = g; 

                // set descent parameter used in Wolfe conditions 
                // SNESLineSearchBTSetAlpha(linesearch, 0.1); 
                
                // linesearch  - The linesearch context
                // X   - The current solution
                // F   - The current function
                // fnorm   - The current norm
                // Y   - The search direction
                SNESLineSearchApply(linesearch, raw_type(x_p), raw_type(g_p), &fnorm, raw_type(d));

                SNESLineSearchGetLambda(linesearch, &lambda); 
                alpha = lambda; 

            }
            else
            {
                Vector x_0 = x, x_k = x;
                Scalar alpha_c, alpha_p, dg = dot(d,g);
                Scalar f, f0, fc, fp, t1, t2, t3, a, b, disc; 
                alpha = 1; 

                if(dg >= 0)
                {
                    if(mpi_rank == 0) {
                        std::cerr<< "utopia::LS::backtracking:: d is not descent direction \n"; 
                        assert(false);
                    }

                    alpha = 0.0;
                    return false;
                }

                fun.value(x_k, f);
                f0 = f;
                fc = f;
                alpha_c = alpha; 

                Scalar it = 0; 

                if(verbose_)
                    PrintInfo::print_init("BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| alpha ||"}); 

                while(alpha > c2_ && it < max_it_)
                {
                    x_k = x_0 + alpha * d;
                    fun.value(x_k, f);

                    // check decrease condition (wolfe condition)
                    if(f < f0 + c1_ * alpha * dg )
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
                    if(verbose_)
                        PrintInfo::print_iter_status({it, alpha}); 
                }

            }

            return true; 
        }

    public: 
        bool set_parameters(const Parameters params) override
        {
            verbose_    = params.line_search_inner_verbose(); 
            c1_         = params.c1(); 
            c2_         = params.c2(); 
            max_it_     = params.n_line_search_iters(); 
            rho_        = params.ls_rho();
            alpha_min_   = params.alpha_min(); 
            
            return true; 
        } 

    private:
        SizeType mpi_size = mpi_world_size();
        SizeType mpi_rank = mpi_world_rank();

        bool verbose_;      /*!< Verbose inside of LS strategy.  */  
        Scalar c1_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */  
        Scalar c2_;         /*!< Constant for Wolfe conditions \f$ c_1 \in (0,1),   c_1 = 10^{-4} \f$.  */  
        Scalar max_it_;     /*!< Maximum of the iterations inside of LS strategy.  */  
        Scalar rho_;        /*!< Contraction factor.   */  
        Scalar alpha_min_;  /*!< Minimum allowed step-size.   */  


    }; 
#endif //WITH_PETSC 

}


#endif //UTOPIA_QUAD_CUB_BACKTRACKING_HPP

