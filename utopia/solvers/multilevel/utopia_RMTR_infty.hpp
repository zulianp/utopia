/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-30
*/

#ifndef UTOPIA_RMTR_INFTY_HPP
#define UTOPIA_RMTR_INFTY_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"


#include "utopia_Linear.hpp"

#include "petscmat.h"
#include "petscvec.h"
#include <petsc/private/snesimpl.h>
#include "petscsnes.h"  

#include "utopia_Level.hpp"


#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"


namespace utopia 
{
    /**
     * @brief      The class for Nonlinear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType>
    class RMTR_infty : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::Transfer<Matrix, Vector>   Transfer;


        typedef utopia::Level<Matrix, Vector>               Level;

    

    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR_infty(    const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
                const std::shared_ptr<Solver> &coarse_solver = std::shared_ptr<Solver>(),
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _smoother(smoother), 
                _coarse_solver(coarse_solver) 
        {
            set_parameters(params); 
        }

        virtual ~RMTR_infty(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            
            _parameters = params; 
        }

        virtual bool solve(FunctionType & fine_fun, Vector &x_h)
        {
            Vector rhs = local_zeros(local_size(x_h)); 
            return solve(fine_fun,  x_h, rhs); 
        }

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(FunctionType &fine_fun, Vector & x_h, const Vector & rhs) 
        {
            this->init_solver("RMTR_infty", {" it. ", "|| r_N ||", "r_norm" , "E"}); 



            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;
            std::cout<<"RMTR_infty: number of levels: "<< l << "  \n"; 

            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 


            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    multiplicative_cycle(fine_fun, x_h, rhs, l); 
                // else if(this->cycle_type() =="full")
                //     full_cycle(rhs, l, x_0); 
                else
                    std::cout<<"ERROR::UTOPIA_MG<< unknown MG type... \n"; 


                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_h) == 1)
                    {
                        x_h = local_zeros(local_size(x_h));
                        return true; 
                    }
                #endif    

                fine_fun.gradient(x_h, F_h); 
                
                Scalar energy; 
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(F_h);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy}); 

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                it++; 
            
            }


            Vector error = hessian * x_h - F_h; 
            Scalar e_norm = norm2(error); 

            std::cout<<"e_norm:  "<< e_norm << " \n"; 

            return true; 
        }



    private: 

        inline FunctionType &levels(const SizeType &l)
        {
            return this->_nonlinear_levels[l]; 
        }

        inline Transfer &transfers(const SizeType & l)
        {
            return this->_transfers[l]; 
        }

        inline Level &lin_levels(const SizeType &l)
        {
            return this->_levels[l]; 
        }



        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {
            Vector g_fine, g_coarse, g_diff, r_h,  r_2h, u_2l, s_coarse, s_fine, u_init; 

           this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, g_fine);             

            r_h = g_fine - f; 


            transfers(l-2).restrict(r_h, r_2h); 
            transfers(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(l-2), u_2l); 

            
            levels(l-2).gradient(u_2l, g_coarse); 


            u_init = u_2l; 
            g_diff = r_2h - g_coarse;  // tau correction 
                          

            if(l == 2)
            {
                coarse_solve(levels(0), u_2l, g_diff, s_coarse); 

            }
            else
            {
                // // recursive call into FAS - needs to be checked 
                // for(SizeType k = 0; k < this->mg_type(); k++)
                // {   
                //     SizeType l_new = l - 1; 
                //     multiplicative_cycle(levels(l-2), u_2l, L_2l, l_new); 
                // }
            }

            // s_coarse = u_2l - u_init; 

           // this->zero_boundary_correction(levels(0), s_coarse); 

            transfers(l-2).interpolate(s_coarse, s_fine);

            //this->zero_boundary_correction(fine_fun, s_fine); 
            
            u_l += s_fine; 

            // POST-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->post_smoothing_steps()); 

            return true; 
        }


        /**
         * @brief      The function invokes smoothing. 
         *
         * @param[in]  A     The stifness matrix. 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool smoothing(Function<Matrix, Vector> &fun,  Vector &x, const Vector &f, const SizeType & nu = 1)
        {
            _smoother->sweeps(nu); 
            _smoother->nonlinear_smooth(fun, x, f); 
            return true; 
        }


        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & g_diff, Vector & s)
        {   
                //_coarse_solver->solve(fun, x, rhs); 
            
                Vector g; 
                fun.gradient(x, g);


                Vector g_0 = g; 
                g +=g_diff; 


                s = 0*x;

                Scalar g_norm = norm2(g); 

                Scalar energy, energy2, energy3; 

                fun.value(x, energy); 
            

                energy3 = energy; // - dot(g_diff, s); 
                energy2 = energy - dot(g_diff, x); 


                // STATISTICS 
                // std::cout<<"g_norm:  "<< g_norm <<  "     linear residual norm: "<< r_norm<< "    E: "<< energy   <<"    E2: "<< energy2   << " \n"; 
                std::cout<<"g_norm:  "<< g_norm <<  "   E3: "<< energy3 << "    E: "<< energy   <<"    E2: "<< energy2   << " \n"; 


            for(auto i = 0; i < 1; i ++)
            {
                auto linear_solver  = std::make_shared< LUDecomposition<Matrix, Vector> >();

            
                Matrix A; 
                fun.hessian(x, A); 
                

                linear_solver->solve(A, -1*g, s);

                Vector r = -1*g - A * s; 
                Scalar r_norm = norm2(r); 


                x +=s;   


                Vector g; 
                fun.gradient(x, g);
                g +=g_diff; 



                g_norm = norm2(g); 




               this->make_iterate_feasible(fun, x); 


                fun.value(x, energy); 
        
                energy2 = energy - dot(g_diff, x); 
                energy3 = energy - dot(g_diff, s); 


                // STATISTICS 
                // std::cout<<"g_norm:  "<< g_norm <<  "     linear residual norm: "<< r_norm<< "    E: "<< energy   <<"    E2: "<< energy2   << " \n"; 
                std::cout<<"g_norm:  "<< g_norm <<  "   E3: "<< energy3 << "    E: "<< energy   <<"    E2: "<< energy2   << " \n"; 
            

            }

            return true; 
        }


    public:
        /**
         * @brief      Function changes direct solver needed for coarse grid solve. 
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return     
         */
        bool change_coarse_solver(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            _coarse_solver = nonlinear_solver; 
            _coarse_solver->set_parameters(_parameters); 
            return true; 
        }



    protected:   
        std::shared_ptr<Smoother>           _smoother;
        std::shared_ptr<Solver>             _coarse_solver;  

    private:
        Parameters                          _parameters; 


    };

}

#endif //UTOPIA_RMTR_INFTY_HPP

