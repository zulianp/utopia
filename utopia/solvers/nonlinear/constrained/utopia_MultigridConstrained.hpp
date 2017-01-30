/*
* @Author: alenakopanicakova
* @Date:   2016-03-29
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-12
*/

#ifndef UTOPIA_MULTIGRID_SUBDIFFERENTIALS_HPP
#define UTOPIA_MULTIGRID_SUBDIFFERENTIALS_HPP
#include "utopia_Smoother.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_Utils.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_ConvergenceReason.hpp"
#include <ctime>


// ----------------------------------------------------WORK IN PROGRESS-----------------------------------------
// needs to be finished, tested and put as part of TR subproblem

namespace utopia 
{
    /**
     * @brief      The class for Linear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class MultigridConstrained : public Multigrid<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector>        Solver;
        typedef utopia::Smoother<Matrix, Vector>            Smoother;
        typedef utopia::MultiLevelBase<Matrix, Vector>      MultiLevelBase;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;

        


    public:

       /**
        * @brief      Multigrid class for constrained problems. 
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        MultigridConstrained(   const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
                                const std::shared_ptr<Solver> &direct_solver = std::shared_ptr<Solver>()) :
                                
                                Multigrid<Matrix, Vector>(smoother, direct_solver)                                                             
        {

        }

        ~MultigridConstrained(){} 
        
        inline Level &levels(const SizeType &l)
        {
            return MultiLevelBase::_levels[l]; 
        }

        inline Transfer &transfers(const SizeType &l)
        {
            return MultiLevelBase::_transfers[l]; 
        }

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param[in]  ub    The upper bound.
         * @param[in]  lb    The lower bound.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(const Vector &rhs, Vector & x_0, const Vector &ub, const Vector &lb)
        {
            SizeType l = this->num_levels(); 
            Vector b = rhs; 

            Scalar r_norm, r0_norm, rel_norm, e_norm, contact_size;
            SizeType it = 0; 

            bool converged = false; 
            
            this->init_solver("Multigrid for constrained problems", {" it. ", "|| r_N ||", "r_norm", "||e||", "# new AP" }); 
            Vector r_h, r_H, c_H, c_h; 

            while(!converged)
            {            

                std::vector<SizeType> active_nodes;    // this requires check
                nonlinear_smoothing(levels(l-1).A(), b, ub, lb, x_0, active_nodes, this->pre_smoothing_steps()); 
                transfers(l-2).apply_truncated_basis_to_interpolation(active_nodes); 
                levels(l-1).enforce_active_set(active_nodes, x_0, b); 

                // finish this properly ... 
                this->galerkin_assembly(levels(l-1).A()); 

                // residual transfer 
                r_h = b - levels(l-1).A() * x_0; 

                if(it == 0)
                    r0_norm = norm2(r_h); 

                transfers(l-2).restrict(r_h, r_H); 

                // prepare correction 
                c_H = zeros(r_H.size().get(0));        


                // call into MG level hierarchy
                this->mg(r_H, l - 1, c_H); 


                // correction transfer
                transfers(l-2).interpolate(c_H, c_h); 
                x_0 += c_h; 


                // -------------------------- statistics --------------------------//
                r_h = b - levels(l-1).A() * x_0; 
                r_norm = norm2(r_h);
                rel_norm = r_norm/r0_norm; 

                e_norm = norm2(c_h); 
                contact_size = active_nodes.size(); 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, e_norm, contact_size}); 

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                it++; 
            
            }
            return true; 
        }


    private:
        /**
         * @brief      The function invokes smoothing with projecting constraints. 
         *
         * @param[in]  A     The stifness matrix. 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool nonlinear_smoothing(const Matrix &A, const Vector &rhs, const Vector& ub, const Vector& lb, Vector &x, std::vector<SizeType>& zero_rows, const SizeType & nu = 1)
        {
            this->_smoother->sweeps(nu); 
            this->_smoother->nonlinear_smooth(A, rhs, ub, lb, x, zero_rows);
            return true; 
        }


    };

}

#endif //UTOPIA_MULTIGRID_SUBDIFFERENTIALS_HPP

