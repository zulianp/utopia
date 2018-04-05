/*
* @Author: kopanicakova
* @Date:   2018-02-02 10:21:12
* @Last Modified by:   kopanicakova
* @Last Modified time: 2018-02-06
*/
#ifndef UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
#define UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP

#include "utopia_Core.hpp" 
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>
#include <cmath>


namespace utopia
{
    enum IterationStatus  { VERY_SUCCESSFUL = 2, 
                            SUCCESSFUL      = 1, 
                            UN_SUCCESSFUL   = 0};


    enum IterationType  {   MONOLITHIC  = 0, 
                            STAG_1      = 1, 
                            STAG_2      = 2};



    template<class Matrix, class Vector, class FunctionType>
    class AutoAdaptiveAlternateMin :    public TrustRegionBase<Matrix, Vector>, 
                                        public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     NonlinearSolver;
        typedef utopia::TRSubproblem<Matrix, Vector>        TRSubproblem; 
        typedef utopia::TrustRegionBase<Matrix, Vector>     TrustRegionBase; 

    public:
       AutoAdaptiveAlternateMin(    const std::shared_ptr<TRSubproblem> &tr_subproblem_stag1 = std::shared_ptr<TRSubproblem>(),
                                    const std::shared_ptr<TRSubproblem> &tr_subproblem_stag2 = std::shared_ptr<TRSubproblem>(),
                                    const std::shared_ptr<TRSubproblem> &tr_subproblem_mono = std::shared_ptr<TRSubproblem>(),
                                    const Parameters params                                 = Parameters() ) :
                                    
                                    NonlinearSolver(tr_subproblem_mono, params), 
                                    _stag1_subproblem(tr_subproblem_stag1), 
                                    _stag2_subproblem(tr_subproblem_stag2),
                                    _monolithic_subproblem(tr_subproblem_mono)
                                    
        {
            set_parameters(params); 
        }



        bool solve( AlternateMinLevel<Matrix, Vector,FunctionType > level, 
                    Vector &x_mono, 
                    Vector &x_stag1, 
                    Vector &x_stag2) 
        {
            using namespace utopia;

            auto stag1_nl_solver = std::make_shared<TrustRegion<Matrix, Vector> >(_stag1_subproblem);
            auto stag2_nl_solver = std::make_shared<TrustRegion<Matrix, Vector> >(_stag2_subproblem);

            bool converged = false;     

            NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

            SizeType global_it = 0, stag1_it = 0; 
            Scalar ared, pred, rho, E_mono_old, E_mono_new; 
            Scalar g_norm, g0_norm, r_norm, s_norm; 
            Scalar global_delta = this->delta0(); 

            IterationStatus iteration_status = VERY_SUCCESSFUL; 
            IterationType   iteration_type   = MONOLITHIC; 
            Vector x_trial_mono = x_mono; 

            Vector g_mono = local_zeros(local_size(x_mono)); 
            Vector s_mono = local_zeros(local_size(x_mono)); 
            level.fun_monolithic().gradient(x_mono, g_mono); 

            Scalar it_scaled=0.0; 
            Scalar eta_inex=0.0, rho_inex =0.0; 

            Matrix H_mono; 
            this->init_solver("AutoAdaptiveAlternateMin", {" it. ", "|| g ||", "J_k", "J_{k+1}", "ared", "pred", "rho", "delta_k", "it_satus", "it_type", "eta_l", "rho_l"}); 

            while(!converged)
            {
                level.fun_monolithic().value(x_mono, E_mono_old); 
                level.fun_monolithic().hessian(x_mono, H_mono); 

                //----------------------------------------------------------------------------
                //     new step p_k w.r. ||p_k|| <= delta. - monolithic one 
                //----------------------------------------------------------------------------      
                if(iteration_type == MONOLITHIC)
                {
                    if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
                        tr_subproblem->current_radius(global_delta);  

                    this->linear_solve(H_mono, g_mono, s_mono);
                    this->get_pred(g_mono, H_mono, s_mono, pred); 

                    x_trial_mono = x_mono + s_mono; 

                    stag1_it = 0; 

                    it_scaled++; 
                }
                else if (iteration_type == STAG_1)
                {
                    //  transfer ...
                    level.transfer_from_monolithic_to_stag1(x_mono, x_stag1); 
                    level.transfer_from_monolithic_to_stag2(x_mono, x_stag2); 
                    level.transfer_from_stag2_to_stag1(x_stag2);  // transfer to aux 

                    Vector g_stag1 = local_zeros(local_size(x_stag1)); 
                    Vector s_stag1 = local_zeros(local_size(x_stag1)); 
                    Matrix H_stag1; 

                    level.fun_stag1().hessian(x_stag1, H_stag1); 
                    level.fun_stag1().gradient(x_stag1, g_stag1); 

                    if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
                        tr_subproblem->current_radius(global_delta);  

                    this->linear_solve(H_stag1, g_stag1, s_stag1);

                    // local trial step 
                    Vector x_trial_stag = x_stag1 + s_stag1; 

                    // interpolate to global... 
                    level.transfer_from_stag1_to_monolithic(x_trial_stag, x_trial_mono); 
                    s_mono = x_trial_mono - x_mono; 

                    this->get_pred(g_mono, H_mono, s_mono, pred); 

                    stag1_it++; 
                    it_scaled += 0.666666666; 

                }
                else if (iteration_type == STAG_2)
                {
                    level.transfer_from_monolithic_to_stag1(x_mono, x_stag1); 
                    level.transfer_from_monolithic_to_stag2(x_mono, x_stag2); 
                    level.transfer_from_stag1_to_stag2(x_stag1);  // transfer to aux 


                    Vector g_stag2 = local_zeros(local_size(x_stag2)); 
                    Vector s_stag2 = local_zeros(local_size(x_stag2)); 
                    Matrix H_stag2; 

                    Scalar energy_old_local, energy_new_local; 
                    level.fun_stag2().value(x_stag2, energy_old_local); 

                    level.fun_stag2().hessian(x_stag2, H_stag2); 
                    level.fun_stag2().gradient(x_stag2, g_stag2); 

                    if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
                        tr_subproblem->current_radius(global_delta);  

                    this->linear_solve(H_stag2, g_stag2, s_stag2);

                    // local trial step 
                    Vector x_trial_stag = x_stag2 + s_stag2; 
                    level.fun_stag2().value(x_trial_stag, energy_new_local); 

                    // interpolate to global... 
                    level.transfer_from_stag2_to_monolithic(x_trial_stag, x_trial_mono); 
                    s_mono = x_trial_mono - x_mono; 

                    this->get_pred(g_mono, H_mono, s_mono, pred); 
                    stag1_it = 0; 

                    it_scaled += 0.33333333333; 
                }
                //----------------------------------------------------------------------------
                //     acceptance of trial point 
                //----------------------------------------------------------------------------
                level.fun_monolithic().value(x_trial_mono, E_mono_new);

                // decrease ratio 
                ared = E_mono_old - E_mono_new;             // reduction observed on objective function
                rho = ared/ std::abs(pred);                 // decrease ratio         


                if(rho >= this->rho_tol())
                    x_mono = x_trial_mono; 

                //----------------------------------------------------------------------------
                //    convergence check 
                //----------------------------------------------------------------------------

                Vector Hs = H_mono * s_mono; 
                Hs += g_mono; 

                level.fun_monolithic().gradient(x_mono, g_mono); 
                Scalar g_norm_old = g_norm; 
                g_norm = norm2(g_mono); 
                
                Scalar Hs_norm = norm2(Hs);

                eta_inex = (g_norm - Hs_norm)/ g_norm_old; 
                rho_inex = (g_norm_old - g_norm)/ (g_norm_old - Hs_norm); 

                global_it++; 
                PrintInfo::print_iter_status(global_it, {g_norm, E_mono_old,  E_mono_new, ared, pred,  rho, global_delta, iteration_status, iteration_type, eta_inex, rho_inex}); 
                converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), global_it, g_norm, 9e9, 9e9, global_delta); 

                //----------------------------------------------------------------------------
                //      tr. radius update 
                //----------------------------------------------------------------------------
                if(rho < this->eta1())
                {
                    global_delta *= this->gamma1();
                    iteration_status = UN_SUCCESSFUL; 
                }
                else if(rho > this->eta2())
                {
                    global_delta = std::min(this->gamma2() * global_delta, this->delta_max()); 
                    iteration_status = VERY_SUCCESSFUL; 
                }
                else
                    iteration_status = SUCCESSFUL; 

                // complete alternation 
                if(rho_inex > 0.8 || iteration_type==STAG_1)
                    iteration_type = MONOLITHIC; 
                else if(rho_inex < 0.8 && iteration_type == MONOLITHIC)
                    iteration_type = STAG_1; 
                else 
                    iteration_type = STAG_2;
            }

            // just to have nice pictures on each ts ... 
            level.transfer_from_monolithic_to_stag1(x_mono, x_stag1); 
            level.transfer_from_monolithic_to_stag2(x_mono, x_stag2); 


            // some benchmarking 
            auto data_path = Utopia::instance().get("auto_adaptive_path");
            if(!data_path.empty())
            {
                CSVWriter writer; 
                if (mpi_world_rank() == 0)
                {
                  if(!writer.file_exists(data_path))
                  {
                      writer.open_file(data_path); 
                      writer.write_table_row<std::string>({("global_it"), "it_scaled"}); 
                  }
                  else
                      writer.open_file(data_path); 
                  
                  writer.write_table_row<Scalar>({Scalar(global_it), Scalar(it_scaled)}); 
                  writer.close_file(); 
                }
            }



            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            NonlinearSolver::set_parameters(params);
            TrustRegionBase::set_parameters(params);
        }


        bool solve(Function<Matrix, Vector> &fun, Vector &x_k) override 
        {
            std::cerr<<"-------AutoAdaptiveAlternateMin::classical solve not implemented ---------  \n"; 
            return false; 
        }


    private: 

        std::shared_ptr<TRSubproblem>        _stag1_subproblem;     
        std::shared_ptr<TRSubproblem>        _stag2_subproblem;   
        std::shared_ptr<TRSubproblem>        _monolithic_subproblem;   




    };

}
#endif //UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
