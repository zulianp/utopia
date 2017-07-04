/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-04
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
     * @brief      The class for Full approximation scheme.
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
        * @brief      The class for Full approximation scheme.
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
        

        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        void set_parameters(const Parameters params) override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            
            _parameters = params; 
        }

        /**
         *
         * @return     Name of class
         */
        virtual std::string name_id() override
        {
            return "FAS"; 
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


        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l) override
        {
            Vector g_fine, g_coarse, u_2l, e, u_init; 

           this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, g_fine);             

            g_fine -= f; 

            transfers(l-2).restrict(g_fine, g_fine); 
            transfers(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(l-2), u_2l); 
            this->zero_correction_related_to_equality_constrain(levels(l-2), g_fine); 

            levels(l-2).gradient(u_2l, g_coarse); 

            u_init = u_2l; 
            g_coarse -= g_fine;  // tau correction 
                                 
            if(l == 2)
            {
                this->make_iterate_feasible(levels(0), u_2l); 
                coarse_solve(levels(0), u_2l, g_coarse); 
            }
            else
            {
                // recursive call into FAS - needs to be checked 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    this->multiplicative_cycle(levels(l-2), u_2l, g_coarse, l_new); 
                }
            }

            e = u_2l - u_init; 
            transfers(l-2).interpolate(e, e);
            this->zero_correction_related_to_equality_constrain(fine_fun, e); 
            
            u_l += e; 

            // POST-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->post_smoothing_steps()); 

            return true; 
        }


        /**
         * @brief      The function invokes smoothing. 
         *
         * @param      fun   Function to be minimized 
         * @param[in]  f     The right hand side.
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



        /**
         * @brief      The function invokes coarse solver. 
         *
         * @param      fun   Function to be minimized 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & rhs) override
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

