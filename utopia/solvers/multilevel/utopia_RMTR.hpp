/*
* @Author: alenakopanicakova
* @Date:   2017-05-04
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-06
*/

#ifndef UTOPIA_RMTR_HPP
#define UTOPIA_RMTR_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_TR_base.hpp"
#include "utopia_TrustRegion.hpp"



// #include <petscksp.h>
// #include <petscsys.h>

// #include "petscmat.h"
// #include "petscvec.h"


namespace utopia 
{
    /**
     * @brief      The class for Line-search multilevel optimization algorithm. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType>
    class RMTR : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>, 
                 public TrustRegionBase<Matrix, Vector> 
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::TrustRegion<Matrix, Vector>         TRSolver;
        typedef utopia::TRSubproblem<Matrix, Vector>            TRSubproblem; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;

    

    public:

        RMTR(const std::shared_ptr<TRSubproblem> &tr_subproblem = std::shared_ptr<TRSubproblem>(), 
            const Parameters params = Parameters()):
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _eps_smooth(0.5e-9), 
                _Kg(0.5), 
                _eps_delta(0.001), 
                _delta0(100),
                _tr_subproblem(tr_subproblem)
                
        {
            set_parameters(params); 
        }

        virtual ~RMTR(){} 
        

        void set_parameters(const Parameters params)  override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            TrustRegionBase<Matrix, Vector>::set_parameters(params);

            // add additional things 
            _parameters = params; 
        }

        virtual bool solve(FunctionType & fine_fun, Vector &x_h)
        {
            std::cout<<"RMTR::solve... \n"; 
            Vector rhs = local_zeros(local_size(x_h)); 
            return this->solve(fine_fun,  x_h, rhs); 
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
            this->init_solver("RMTR", {" it. ", "|| r_N ||", "r_norm" , "E"}); 



            fine_fun.get_boundary_values(x_h); 


            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;

            std::cout<<"RMTR: number of levels: "<< l << "  \n"; 

            // just to check what is problem 
            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 

            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            while(!converged)
            {            
                //if(this->cycle_type() =="multiplicative")
                    multiplicative_cycle(fine_fun, x_h, rhs, l); 


                // else if(this->cycle_type() =="full")
                //     full_cycle(rhs, l, x_0); 
                // else
                //     std::cout<<"ERROR::MG_OPT<< unknown MG type... \n"; 


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
                // TODO:: change 
                converged = NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::check_convergence(it, r_norm, rel_norm, 1); 
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
            std::cout<<"-------- at some point: RMTR ....  multiplicative_cycle ....... \n"; 

            Scalar delta = _delta0, delta_local = _delta0; 
            Scalar  intermediate_delta = delta; 

            Scalar rho, ared, pred; 

            Vector g_l, g_2l, u_2l, e_2h, e_h, u_init; 

            Matrix H_l, H_2l; 

            Scalar E_k_h, E_k1_h;       // energy on fine level 
            Scalar E_k_2h, E_k1_2h;     // energy on coarse level 




        for(int it = 0; it < 5000; it ++)
        {

            fine_fun.gradient(u_l, g_l);  

            Vector smooth_s =  0 * u_l; 
            Vector u_l_init = u_l; 
                 
            delta = 10000; 

            local_smoother_TR(fine_fun, u_l, g_l, delta, 0, l, pred, smooth_s); 


            fine_fun.gradient(u_l, g_l);  
            Scalar r_norm = norm2(g_l); 
            bool converged = NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::check_convergence(it, r_norm, 1, 1); 

            if(converged)
                exit(1); 


            fine_fun.value(u_l, E_k_h);

            transfers(l-2).restrict(g_l, g_2l);    
            transfers(l-2).project_down(u_l, u_2l);   

            Scalar delta_fine = delta; 

            Vector u_2l_init = u_2l; 
               
            local_smoother_TR(levels(0), u_2l, g_2l, delta, 1, l-1, pred, e_2h); 

            std::cout<<"dot(g_2l, e_2h): "<<dot(g_2l, e_2h) << "    \n"; 


            // here we take initial delta from fine level 
            delta = delta_fine; 

            // correction is up
            transfers(l-2).interpolate(e_2h, e_h);
            // this->zero_boundary_correction(fine_fun, e_h); 


            std::cout<<"dot(g_l, e_h): "<<dot(g_l, e_h) << "    \n"; 


            Vector u_trial = u_l + e_h; 
            fine_fun.value(u_trial, E_k1_h);
        
            ared = E_k_h - E_k1_h;  


            if(ared <= 0 || pred <= 0)
            {
                rho = 0; 
            }
            else if(ared == pred)
            {
                rho = 1; 
            } 
            else
            {
                // agrement between glob and local function
                rho =  ared/pred; 
            }

            std::cout<<"-----------RMTR---------    rho: "<< rho << "    ared:  "<< ared <<  "    pred: "<< pred << "  \n"; 


            //----------------------------------------------------------------------------
            //     acceptance of trial point 
            //----------------------------------------------------------------------------


             // good reduction, accept trial point 
              if (rho >= this->rho_tol())
              {
                std::cout<<"RMTR:: TAKING point from previous level ... \n"; 
                u_l = u_trial; 

              }
              // otherwise, keep old point
              else
              {
                std::cout<<"RMTR:: NOT  taking point from previous level ... \n"; 
              }


            // _______________________________ NO RECURSION YET _______________________________

            //----------------------------------------------------------------------------
            //      tr. radius update 
            //----------------------------------------------------------------------------
            

                if(rho < this->eta1())
                {   
                     delta = this->gamma1() * delta; 
                }
                else if (rho > this->eta2() )
                {
                     delta = this->gamma2() * delta; 
                }
                else
                {
                    delta = delta; 
                }     

                if(delta < 1e-12)
                    exit(1); 




            }

            return true; 
        }











    public:


    protected:   




    private:


        Scalar level_dependent_norm(const Vector & u, const SizeType & current_l)
        {
            // if(current_l == this->num_levels())
            // {
            //     return 0.0; 
            // }
            // else
            // {
            //     Vector s = u; // carries over prolongated correction
         
            //     for(SizeType i = current_l; i < this->num_levels(); i++)
            //     {
            //         // std::cout<<"yes, sum ... i:  "<< i << " \n"; 
            //         transfers(i-1).interpolate(s, s); 
            //     }
            //     //std::cout<<"level_dependent_norm: "<< norm2(s) << "   \n"; 
            //     return norm2(s); 
            // }
            

            if(current_l == this->num_levels())
            {
                return 0.0; 
            }
            else
            {
                // Vector s = u; // carries over prolongated correction
                Vector s = u; // carries over prolongated correction
         
                for(SizeType i = current_l; i < this->num_levels(); i++)
                {
                    transfers(i-1).interpolate(s, s); 
                    transfers(i-1).restrict(s, s); 
                    // s = dot(s,s); 
                }
             //   std::cout<<"level_dependent_norm graffon :   "<< norm2(s) << "   \n"; 
                return norm2(s); 
            }

        }



        bool local_smoother_TR(FunctionType & fun, Vector & x_k, const Vector & Rg, Scalar & delta, const SizeType & max_it, const SizeType & level, Scalar & pred_local_model, Vector & s_global)
        {
            const Scalar delta0 = delta; 
            this->make_iterate_feasible(fun, x_k); 


            // check on smoothness of restrcited gradient 
            Scalar smoothness_norm = norm2(Rg); 
            if(smoothness_norm < _eps_smooth && level<this->num_levels())
            {
                std::cout<<"local_smoother_TR:: termination since beginning based on smoothness.... \n"; 
                pred_local_model = 0; 
                std::cout<<"pred: "<< pred_local_model << " \n"; 
                return true; 
            }

            Scalar s_norm, g_norm, intermediate_delta = delta0; 
            Scalar E, E_k1, E_k; 
            Scalar ared, pred, rho; 

            Vector g, g_delta, s = x_k, x_init = x_k; 
            Matrix H, H_delta; 


            SizeType it_success = 0, it = 0;
            bool terminate = false;  


            // first order consistency 
            fun.gradient(x_k, g);
            g_delta = Rg - g;   


            std::cout<<"||g_delta|| "<< norm2(g_delta) << " \n"; 


            // in order to measure model decrease 
            g_norm = norm2(g + g_delta); 
            fun.value(x_k, E); 
            Scalar E0 = E + dot(g_delta, x_k); // + 0.5 * dot(x_k, H_delta * x_k); 


            if(level == 1)
            {
                this->init_solver("COARSE SOLVE", {" it. ", "|| g ||", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred", 
                                            "rho", "delta_k", "|| p_k || "});
            }
            else
            {
                this->init_solver("SMOOTHER", {" it. ", "|| g ||", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred", 
                                            "rho", "delta_k", "|| p_k || "});
            }


            PrintInfo::print_iter_status(0, {g_norm, 0, 0, 0, 0, 0, 0, delta, 0}); 



            while(!terminate)
            {
                
                fun.gradient(x_k, g);
                g = g + g_delta; 

                fun.value(x_k, E_k); 
                E_k = E_k + dot(g_delta, x_k); 

                fun.hessian(x_k, H); 


                if(level >1)
                {
                    Vector sg = -1 * g; 
                    Vector s_smooth; 

                    std::cout<<"smoothing.... \n"; 
                    
                    // MatSOR( raw_type(H), 
                    // raw_type(sg), 
                    // 1, 
                    // //  SOR_FORWARD_SWEEP,
                    // SOR_LOCAL_SYMMETRIC_SWEEP,    // parallel implementation - builds block jacobi and on blocks it calls GS 
                    // 0, 
                    // 2, 
                    // 2, 
                    // raw_type(s)); 

                    this->get_pred(g, H, s, pred); 
                    std::cout<<"norm2(s):   " << norm2(s) << "  \n"; 

                }
                else
                {

                    // strange thing 
                    // KSP ksp; 
                    // PC pc; 
                    // MPI_Comm            comm; 
                    // PetscObjectGetComm((PetscObject)raw_type(H), &comm);
                    // KSPCreate(comm, &ksp);

                    // KSPSetOperators(ksp, raw_type(H), raw_type(H));
                    // KSPSetType(ksp, KSPSTCG);  
                    // KSPGetPC(ksp, &pc);
                    // PCSetType(pc, PCASM);
                    // KSPSetTolerances(ksp,1e-15, 1e-15, PETSC_DEFAULT, PETSC_DEFAULT); 
                    // KSPSTCGSetRadius(ksp, 10000000);
                    // KSPSolve(ksp, raw_type(g),raw_type(s));
                    // KSPSTCGGetObjFcn(ksp, &pred);

                    // since Newton iteration is defined with - 
                    pred = -pred; 
                    s *=-1;  
                }

            //----------------------------------------------------------------------------
            //     building trial point 
            //----------------------------------------------------------------------------
                
                Vector tp = x_k + s;  
                fun.value(tp, E_k1);
                E_k1 = E_k1 + dot(g_delta, tp); // + 0.5 * dot(tp, H_delta * tp); ; 

                // decrease ratio 
                ared = E_k - E_k1;                  // reduction observed on objective function
                pred = std::abs(pred);     

            //----------------------------------------------------------------------------
            //     acceptance of trial point 
            //----------------------------------------------------------------------------
              
            if(ared <= 0 || pred <= 0)
            {
                rho = 0; 
            }
            else if(ared == pred)
            {
                rho = 1; 
            }
            else
            {
               // decrease ratio     
                rho =  ared/pred; 
            } 

              // good reduction, accept trial point 
              if (rho >= this->rho_tol())
              {
                x_k = tp; 
                E = E_k1; 
                it_success++; 
              }
              // otherwise, keep old point
              else
              {
                E = E_k; 
              }

            //----------------------------------------------------------------------------
            //    convergence check 
            //----------------------------------------------------------------------------
              fun.gradient(x_k, g); 
              g = g + g_delta; 
              g_norm = norm2(g); 
              s_norm = norm2(s); 

              if(it_success>= max_it || g_norm <= 1e-7)
              {
                terminate = true; 
              }

            //----------------------------------------------------------------------------
            //    CHECK FOR SMOOTHNESS OF GIVEN LEVEL
            //----------------------------------------------------------------------------

            // compute global correction
            // this could be done by summing up corrections 
            s_global = x_k - x_init; 


            Scalar smoothness_norm = norm_infty(g); 
            if(smoothness_norm < _eps_smooth && level<this->num_levels())
            {
                PrintInfo::print_iter_status(it, {g_norm, E, E_k, E_k1, ared, pred, rho, delta, s_norm}); 
                std::cout<<"termination based on smoothness of grad .... \n"; 
                pred_local_model = E0 - E; 
                std::cout<<"pred: "<< pred_local_model << " \n"; 
                return true; 
            }


            // we need to guarantee that iterates stay in TR radius from previous level 
            Scalar corr_norm = level_dependent_norm(s_global, level);
            if(corr_norm > (1 - _eps_delta)* delta0 && level<this->num_levels())
            {
                PrintInfo::print_iter_status(it, {g_norm, E, E_k, E_k1, ared, pred, rho, delta, s_norm}); 
                std::cout<<"termination due to too long correction .... \n"; 
                pred_local_model = E0 - E; 
                return true; 
            }


            //----------------------------------------------------------------------------
            //      tr. radius update 
            //----------------------------------------------------------------------------
                if(rho < this->eta1())
                {   
                     intermediate_delta = this->gamma1() * delta; 
                }
                else if (rho > this->eta2() )
                {
                     intermediate_delta = this->gamma2() * delta; 
                }
                else
                {
                    intermediate_delta = delta; 
                }      

                // on the finest level we work just with one radius 
                if(level==this->num_levels())
                {
                    delta = intermediate_delta; 
                }
                else
                {
                    // compute radius 
                    Scalar delta_help = delta0 - corr_norm; 
                    delta = std::min(intermediate_delta, delta_help); 
                }


                if(delta < 1e-12)
                {   std::cout<<"EROROR:: delta smaller than 0 ... something is wrong ?? \n ";              
                    pred_local_model = E0 - E; 
                    return true; 
                }

                it++; 
                PrintInfo::print_iter_status(it, {g_norm, E, E_k, E_k1, ared, pred, rho, delta, s_norm}); 
            }

            pred_local_model = E0 - E; 
            std::cout<<"pred: "<< pred << " \n"; 

            return true; 
        }






    private:
        Parameters                          _parameters; 

        Scalar _eps_smooth; 
        Scalar _Kg; 

        Scalar _eps_delta; 
        Scalar _delta0;                     // initial tr radius - restr from fine level 


        std::shared_ptr<TRSubproblem> _tr_subproblem;  // solving TR subproblems on different levels   


    };

}

#endif //UTOPIA_RMTR_HPP

