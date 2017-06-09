/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-09
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

    

        enum Coherence_level {  FIRST_ORDER  = 1, 
                                SECOND_ORDER = 2, 
                                GALERKIN     = 0};


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
                _coarse_solver(coarse_solver), 
                _coherence(FIRST_ORDER) 
        {
            set_parameters(params); 
        }

        virtual ~RMTR_infty(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            
            _delta_init = 1; 

            // _delta.resize(this->num_levels()); 

            // for(Scalar i = 0; i < this->num_levels(); i ++)
            //     _delta.push_back(_delta_init); 

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
            this->init_solver("RMTR_infty", {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 



            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;
            std::cout<<"RMTR_infty: number of levels: "<< l << "  \n"; 

            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 


            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            this->make_iterate_feasible(fine_fun, x_h); 
            Scalar energy; 
            fine_fun.value(x_h, energy); 

            Scalar ared=0, pred=0, rho=0; 

            PrintInfo::print_iter_status(it, {r0_norm, energy, ared, pred, rho, 0.0}); 
            it++; 

            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    multiplicative_cycle(fine_fun, x_h, rhs, l, ared, pred, rho); 
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
                
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(F_h);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    // PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy}); 
                    PrintInfo::print_iter_status(it, {r_norm, energy, ared, pred, rho, 0.0}); 

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



        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l, Scalar & ared, Scalar &coarse_reduction, Scalar &rho)
        {
            Vector g_fine, g_coarse, g_diff, r_h,  g_restricted, u_2l, s_coarse, s_fine, u_init; 


            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, g_fine);   



            r_h = g_fine - f; 

            transfers(l-2).restrict(r_h, g_restricted);
            transfers(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(0), u_2l); 

            // here, grad should be zero by default on places where are BC conditions             
            levels(l-2).gradient(u_2l, g_coarse); 

            u_init = u_2l; 

            if(_coherence != GALERKIN)
                this->zero_boundary_correction(levels(0), g_restricted); 


            g_diff = g_restricted - g_coarse;  // tau correction 


            Matrix H_fine, H_restricted, H_coarse, H_diff; 

            if(_coherence == SECOND_ORDER || _coherence == GALERKIN)
            {
                fine_fun.hessian(u_l, H_fine);   
                transfers(l-2).restrict(H_fine, H_restricted);
                
                if(_coherence == SECOND_ORDER)
                    this->zero_boundary_correction_mat(levels(0), H_restricted); 

                levels(l-2).hessian(u_2l, H_coarse); 
                H_diff = H_restricted - H_coarse; 
            }


            if(l == 2)
            {
                coarse_solve(levels(0), u_2l, g_diff, H_diff, s_coarse, coarse_reduction); 
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


            transfers(l-2).interpolate(s_coarse, s_fine);
            this->zero_boundary_correction(fine_fun, s_fine); 

            Scalar E_old, E_new; 
            fine_fun.value(u_l, E_old); 

            Vector u_t = u_l; 

            u_t += s_fine; 
            fine_fun.value(u_t, E_new); 

            ared = E_old - E_new; 
            rho = ared / coarse_reduction; 

           std::cout<<"E_old:  "<< E_old << "  E_new:  "<< E_new << "  \n"; 

            // TODO:: deltas, line-search, 

            if(coarse_reduction<=0)
                rho = 0; 

            if(rho > 0)
                u_l = u_t; 
            else{
                // u_l = u_t; 
                std::cout<<"RMTR:: not taken trial point... \n"; 
            }



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
            

                // Scalar delta = 1000; 

                // Vector g, s;

                // s = 0 * x;  
                // fun.gradient(x, g);

                // Matrix A; 
                // fun.hessian(x, A);


                // KSP ksp; 
                // PC pc; 
                // MPI_Comm            comm; 
                // PetscObjectGetComm((PetscObject)raw_type(A), &comm);
                // KSPCreate(comm, &ksp);

                // KSPSetOperators(ksp, raw_type(A), raw_type(A));
                // KSPSetType(ksp, KSPSTCG);  
                // KSPGetPC(ksp, &pc);
                // PCSetType(pc, PCASM);
                // KSPSetTolerances(ksp,1e-15, 1e-15, PETSC_DEFAULT, 10000); 
                // KSPSTCGSetRadius(ksp, delta);
                // KSPSolve(ksp, raw_type(g),raw_type(s));
                // // KSPSTCGGetObjFcn(ksp, &pred);


                // x = x - s; 



            return true; 
        }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & g_diff, const Matrix & H_diff, Vector & s_global, Scalar & reduction)
        {   
            

                if(_coherence == FIRST_ORDER)
                    std::cout<<"-------- Yes, first order....... \n"; 
                else if(_coherence == SECOND_ORDER)
                    std::cout<<"-------- Yes, second order....... \n"; 


                Scalar energy_old, energy, energy_init, g_norm;
                Scalar ared = 0 , pred = 0, rho = 0; 

                ColorModifier red(FG_LIGHT_YELLOW);
                ColorModifier def(FG_DEFAULT);
                std::cout << red; 

                this->init_solver("COARSE SOLVE", {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 

                Vector x_init = x; 
                Scalar delta = 100000; 
                Vector s = 0 * x; 
                s_global = s; 
                reduction = 0.0; 


                Vector g; 
                Matrix A; 

                
                // --------------------------------------- computation of grad -------------------------------
                if(_coherence != GALERKIN)
                    fun.gradient(x, g);
 
                if(_coherence == FIRST_ORDER)
                    g += g_diff; 
                else if(_coherence == SECOND_ORDER)
                    g += g_diff + H_diff * s_global; 
                else if(_coherence == GALERKIN)
                    g = g_diff + H_diff * s_global; 
                // --------------------------------------------------------------------------------------------


                g_norm = norm2(g); 

                // --------------------------------------- computation of energy -------------------------------
                if(_coherence != GALERKIN)
                    fun.value(x, energy_init); 

                if(_coherence == FIRST_ORDER)
                    energy_init += dot(g_diff, s_global); 
                else if(_coherence == SECOND_ORDER)
                    energy_init += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                else if(_coherence == GALERKIN)
                    energy_init = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                // --------------------------------------------------------------------------------------------
                


                PrintInfo::print_iter_status(0, {g_norm, energy_init, ared, pred, rho, delta }); 



            for(auto i = 0; i < 50; i ++)
            {
                

        
                // --------------------------------------- computation of hessian -------------------------------
                if(_coherence != GALERKIN)
                    fun.hessian(x, A); 

                if(_coherence == SECOND_ORDER)
                    A = A + H_diff; 
                else if(_coherence == GALERKIN)
                    A = H_diff; 
                // --------------------------------------------------------------------------------------------

                // --------------------------------------- computation of energy -------------------------------
                if(_coherence != GALERKIN)
                    fun.value(x, energy_old); 

                if(_coherence == FIRST_ORDER)
                    energy_old += dot(g_diff, s_global); 
                else if(_coherence == SECOND_ORDER)
                    energy_old += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                else if(_coherence == GALERKIN)
                    energy_old = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                // --------------------------------------------------------------------------------------------


                s = 0 * x;


                // strange thing 
                KSP ksp; 
                PC pc; 
                MPI_Comm            comm; 
                PetscObjectGetComm((PetscObject)raw_type(A), &comm);
                KSPCreate(comm, &ksp);

                KSPSetOperators(ksp, raw_type(A), raw_type(A));
                KSPSetType(ksp, KSPSTCG);  
                KSPGetPC(ksp, &pc);
                PCSetType(pc, PCASM);
                KSPSetTolerances(ksp,1e-15, 1e-15, PETSC_DEFAULT, 1000000); 
                KSPSTCGSetRadius(ksp, delta);
                KSPSolve(ksp, raw_type(g),raw_type(s));
                KSPSTCGGetObjFcn(ksp, &pred);

                // since Newton iteration is defined with - 
                pred = -pred; 
                s *=-1;  

            //----------------------------------------------------------------------------
            //     building trial point 
            //----------------------------------------------------------------------------
                
                Vector tp = x + s;  
                s_global += s;   

                // fun.value(tp, energy); 
                // // energy += dot(g_diff, s_global); 
                // energy += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                

                // --------------------------------------- computation of energy -------------------------------
                if(_coherence != GALERKIN)
                    fun.value(tp, energy); 

                if(_coherence == FIRST_ORDER)
                    energy += dot(g_diff, s_global); 
                else if(_coherence == SECOND_ORDER)
                    energy += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                else if(_coherence == GALERKIN)
                    energy = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                // --------------------------------------------------------------------------------------------


                ared = energy_old - energy; 

                // choice for the moment 
                rho = ared/pred; 

                //----------------------------------------------------------------------------
                //     acceptance of trial point 
                //----------------------------------------------------------------------------
                  

                  // good reduction, accept trial point 
                  if (rho >= 0.0000001)
                  {
                    x = tp; 
                    reduction += ared; 
//                  it_success++; 
                  }
                  else
                  {
                   // x = tp; 
                    // since point was not taken 
                    // s_global -= s; 
                  }



                if(rho < 0.1)
                {   
                     delta = 0.2 * delta; 
                }
                else if (rho > 0.85 )
                {
                     delta = 2 * delta; 
                }
                else
                {
                    delta = delta; // just to remember it - delete it later
                }      



                // --------------------------------------- computation of grad -------------------------------
                if(_coherence != GALERKIN)
                    fun.gradient(x, g);
 
                if(_coherence == FIRST_ORDER)
                    g += g_diff; 
                else if(_coherence == SECOND_ORDER)
                    g += g_diff + H_diff * s_global; 
                else if(_coherence == GALERKIN)
                    g = g_diff + H_diff * s_global; 
                // --------------------------------------------------------------------------------------------


                g_norm = norm2(g); 

                PrintInfo::print_iter_status(i, {g_norm, energy, ared, pred, rho, delta}); 

                if(g_norm < 1e-8)
                {
                    std::cout<< def; 
                    return 0; 
                }


            }

            // TODO:: check this in nonlinear case 
            // s = x - x_init; 
            

            std::cout<< def; 

//            exit(1);

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


        bool set_delta(const SizeType & level, const Scalar & radius)
        {
            _delta[level] = radius; 
            return true; 
        }

        Scalar get_delta(const SizeType & level) const 
        {
            return _delta[level]; 
        }


    protected:   
        std::shared_ptr<Smoother>           _smoother;
        std::shared_ptr<Solver>             _coarse_solver;  

        Scalar                              _delta_init; 

        std::vector<Scalar>                 _delta;
        Coherence_level                     _coherence; 

    private:
        Parameters                          _parameters; 


    };

}

#endif //UTOPIA_RMTR_INFTY_HPP

