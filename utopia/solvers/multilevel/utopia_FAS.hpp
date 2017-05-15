/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-15
*/

#ifndef UTOPIA_FAS_HPP
#define UTOPIA_FAS_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"



namespace utopia 
{
    /**
     * @brief      The class for Nonlinear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType>
    class FAS : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::Transfer<Matrix, Vector>   Transfer;

    

    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        FAS(    const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
                const std::shared_ptr<Solver> &coarse_solver = std::shared_ptr<Solver>(),
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _smoother(smoother), 
                _coarse_solver(coarse_solver) 
        {
            set_parameters(params); 
        }

        virtual ~FAS(){} 
        

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
            this->init_solver("FAS", {" it. ", "|| r_N ||", "r_norm" , "E"}); 

            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;

            std::cout<<"FAS: number of levels: "<< l << "  \n"; 

            // just to check what is problem 
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




        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {
            Vector L_l, L_2l, r_h,  r_2h, u_2l, e_2h, e_h, u_init; 
            
            this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, L_l);             

            r_h = L_l - f; 


            transfers(l-2).restrict(r_h, r_2h); 
            transfers(l-2).project_down(u_l, u_2l); 

            
            levels(l-2).gradient(u_2l, L_2l); 

            u_init = u_2l; 
            L_2l = L_2l - r_2h;  // tau correction 
                                 
            if(l == 2)
            {
                this->make_iterate_feasible(levels(0), u_2l); 
                coarse_solve(levels(0), u_2l, L_2l); 
            }
            else
            {
                // recursive call into FAS - needs to be checked 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    multiplicative_cycle(levels(l-2), u_2l, L_2l, l_new); 
                }
            }

            e_2h = u_2l - u_init; 
            transfers(l-2).interpolate(e_2h, e_h);

            this->zero_correction_contributions(fine_fun, e_h); 
            
            u_l += e_h; 

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


        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & rhs)
        {   
            _coarse_solver->solve(fun, x, rhs); 
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

#endif //UTOPIA_FAS_HPP

