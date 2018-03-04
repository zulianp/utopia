/*
* @Author: alenakopanicakova
* @Date:   2017-04-24
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-04
*/

#ifndef UTOPIA_NMGM_HPP
#define UTOPIA_NMGM_HPP
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
    class NonLinearMultigrid : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;

    

    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        NonLinearMultigrid( const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
                            const std::shared_ptr<Solver> &coarse_solver = std::shared_ptr<Solver>(),
                            const Parameters params = Parameters()): 
                            NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                            _smoother(smoother), 
                            _coarse_solver(coarse_solver) 
        {
            set_parameters(params); 
        }

        virtual ~NonLinearMultigrid(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            
            _parameters = params;         
            _sigma      = params.sigma(); 
            this->cycle_type(NESTED_ITERATION); 
        }


        virtual std::string name_id()  override
        {
            return "Nonlinear Multigrid"; 
        }


        virtual bool solve(FunctionType & fine_fun, Vector &x_h) override
        {
            Vector rhs = local_zeros(local_size(x_h)); 
            return solve(fine_fun,  x_h, rhs); 
        }

        /**
         * @brief       The solve function for Nonlinear multigrid method. 
         *              This one is slightly different than in original base class, due to storing vectors from previous iteration
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(FunctionType &fine_fun, Vector & x_h, const Vector & rhs) override
        {
            this->init_solver("Nonlinear multigrid", {" it. ", "|| r_N ||", "r_norm" }); 
            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 

            std::cout<<"NMG: number of levels: "<< l << "  \n"; 
            Scalar r_norm, r0_norm, rel_norm;

            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            // we assume NLMG to be part of nested iteration by default
            // as in "Multi-Grid Methods and Applications, W. Hackbush" 
            std::vector<Vector> rhss; 
            std::vector<Vector> initial_iterates; 
            nested_iteration_cycle(x_h, rhs, l, rhss, initial_iterates); 

            while(!converged)
            {            
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE)
                    multiplicative_cycle(fine_fun, x_h, rhs, l); 
                else if(this->cycle_type() ==NESTED_ITERATION)
                    NMGM(fine_fun, x_h, rhs, l, rhss, initial_iterates); 

                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_h) == 1)
                    {
                        x_h = local_zeros(local_size(x_h));
                        return true; 
                    }
                #endif    

                fine_fun.gradient(x_h, F_h); 
                r_norm = norm2(F_h);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm}); 

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



        // TODO:: make this nicer ! 
        bool nested_iteration_cycle(Vector & u_l, const Vector &/*f*/, const SizeType & l, std::vector<Vector> & rhss, std::vector<Vector> & initial_iterates)
        {
            for(SizeType i = l-2; i >=0; i--)
            {
              transfers(i).restrict(u_l, u_l); 
              this->make_iterate_feasible(levels(i), u_l); 
            }

            Vector L_l = local_zeros(local_size(u_l));
            coarse_solve(levels(0), u_l, L_l); 

            initial_iterates.push_back(u_l); 
            levels(0).gradient(u_l, L_l); 
            rhss.push_back(L_l); 
            
            transfers(0).interpolate(u_l, u_l); 

            for(SizeType i = 1; i <l-1; i++)
            {
                for(SizeType j = 0; j < this->v_cycle_repetition(); j++)
                {
                  Vector f = local_zeros(local_size(u_l));
                  NMGM(levels(i),  u_l, f, i+1, rhss, initial_iterates); 
            
                  initial_iterates.push_back(u_l); 
                  levels(i).gradient(u_l, L_l); 
                  rhss.push_back(L_l); 
                }
                  transfers(i).interpolate(u_l, u_l); 
            }
            return true; 
        }


        bool NMGM(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l, const std::vector<Vector> & rhss, const std::vector<Vector> & initial_iterates)
        {
            Vector L_l, L_2l, r_h,  r_2h, u_2l, e_2h, e_h; 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, L_l);             

            r_h = L_l - f; 
            transfers(l-2).restrict(r_h, r_2h); 
          
            u_2l    = initial_iterates[l-2]; 
            Scalar s = scaling_factor(r_2h); 

            this->zero_correction_related_to_equality_constrain(levels(l-2), r_2h); 
            L_2l = rhss[l-2] - s *r_2h;  // tau correction 
          
            if(l == 2)
            {
                coarse_solve(levels(0), u_2l, L_2l); 
            }
            else
            {
                // recursive call into FAS - needs to be checked 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    NMGM(levels(l-2), u_2l, L_2l, l_new, rhss, initial_iterates); 
                }
            }

            e_2h = 1/s * (u_2l - initial_iterates[l-2]); 
            transfers(l-2).interpolate(e_2h, e_h);

            this->zero_correction_related_to_equality_constrain(fine_fun, e_h); 
            u_l += e_h; 

            // POST-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->post_smoothing_steps()); 
            return true; 
        }

        /**
         * @brief      We keep also this one, since we would like to use it from full-cycle. 
         *
         * @param      fine_fun  Function to be minimized
         * @param      u_l       Iterate
         * @param[in]  f         Rhs 
         * @param[in]  l         Level
         *
         */
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
            Scalar s = scaling_factor(g_fine); 


            u_init = u_2l; 
            g_coarse -= g_fine;  // tau correction 
                                 
            if(l == 2)
            {
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
            
            u_l += 1/s * e; 

            // POST-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->post_smoothing_steps()); 

            return true; 
        }



        /**
         * @brief      Scaling factor 
         *
         * @param[in]  d    Vector of correction. 
         *
         * @return     scaling factor
         */
        Scalar scaling_factor(const Vector & d)
        {
            return _sigma/norm2(d); 
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

        // with RHS
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


        bool set_sigma(const Scalar & sigma)
        {
            _sigma = sigma; 
            return true; 
        }


        Scalar get_sigma() const
        {
            return _sigma; 
        }



    protected:   
        std::shared_ptr<Smoother>           _smoother;
        std::shared_ptr<Solver>             _coarse_solver;  

    private:
        Parameters                          _parameters; 
        Scalar                              _sigma; 

    };

}

#endif //UTOPIA_NMGM_HPP

