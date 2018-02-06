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
    template<class Matrix, class Vector, class FunctionType>
    class AutoAdaptiveAlternateMin :    /*public AlternateMinimization<Matrix, Vector>*/
                                        public TrustRegionBase<Matrix, Vector>, 
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

            SizeType global_it = 0; 
            Scalar ared, pred, rho, E_mono_old, E_mono_new; 
            Scalar g_norm, g0_norm, r_norm, s_norm; 
            Scalar global_delta = this->delta0(); 
            std::cout<<"global_delta: "<< global_delta << " \n";  



            Vector g_mono = local_zeros(local_size(x_mono)); 
            Vector s_mono = local_zeros(local_size(x_mono)); 
            level.fun_monolithic().gradient(x_mono, g_mono); 

            Matrix H_mono; 
            this->init_solver("AutoAdaptiveAlternateMin", {" it. ", "|| g ||", "J_k", "J_{k+1}", "rho", "delta_k"}); 

            
            while(!converged)
            {

                level.fun_monolithic().value(x_mono, E_mono_old); 
                level.fun_monolithic().hessian(x_mono, H_mono); 

                //----------------------------------------------------------------------------
                //     new step p_k w.r. ||p_k|| <= delta. - monolithic one 
                //----------------------------------------------------------------------------      
                if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
                    tr_subproblem->current_radius(global_delta);  

                this->linear_solve(H_mono, g_mono, s_mono);
                this->get_pred(g_mono, H_mono, s_mono, pred); 


                //----------------------------------------------------------------------------
                //     acceptance of trial point 
                //----------------------------------------------------------------------------

                level.fun_monolithic().value(x_mono + s_mono, E_mono_new);

                // decrease ratio 
                ared = E_mono_old - E_mono_new;             // reduction observed on objective function
                rho = ared/ std::abs(pred);                 // decrease ratio         


                if(rho >= this->rho_tol())
                    x_mono += s_mono; 


                //----------------------------------------------------------------------------
                //    convergence check 
                //----------------------------------------------------------------------------

                level.fun_monolithic().gradient(x_mono, g_mono); 
                g_norm = norm2(g_mono); 
                global_it++; 
                converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), global_it, g_norm, 9e9, 9e9, global_delta); 



                PrintInfo::print_iter_status(global_it, {g_norm, E_mono_old,  E_mono_new, rho, global_delta}); 



                //----------------------------------------------------------------------------
                //      tr. radius update 
                //----------------------------------------------------------------------------
                if(rho < this->eta1())
                    global_delta *= this->gamma1();
                else if(rho > this->eta2())
                    global_delta = std::min(this->gamma2() * global_delta, this->delta_max()); 



            }







            // Vector x_mono_test = x_mono; 
            // this->_nl_solver_master->solve(level.fun_monolithic(), x_mono_test); 




            // // ------------------------ staggered --------------------
            // std::cout<<"------------- staggered ---------------------- \n"; 
            // Scalar energy_0, energy_prev, energy_cur; 

            // level.fun_stag1().value(x_stag1, energy_0); 
            // energy_prev = energy_0; 


            // Scalar beta    = 9e9; 
            // SizeType it     = 1; 

            // while(beta > this->get_energy_slope_tol() && it < this->get_num_alternate_steps() )
            // {
            //     this->_nl_solver_master->max_it(1); 
            //     this->_nl_solver_master->solve(level.fun_stag1(), x_stag1); 
            //     level.transfer_from_stag1_to_stag2(x_stag1); 

            //     this->_nl_solver_slave->solve(level.fun_stag2(), x_stag2); 
            //     level.transfer_from_stag2_to_stag1(x_stag2); 


            //     level.fun_stag1().value(x_stag1, energy_cur); 

            //     Scalar val = ((energy_prev - energy_cur)/ (energy_0 - energy_cur)) * it; 
            //     beta = std::atan (val) * 180.0 / M_PI;

            //     energy_prev = energy_cur; 
            //     it++; 


            //     level.transfer_from_stag1_to_monolithic(x_stag1, x_mono); 
            //     level.transfer_from_stag2_to_monolithic(x_stag2, x_mono); 

            //     // ------------------------ monolithic  --------------------
            //     std::cout<<"------------- monolithic ---------------------- \n"; 
            //     x_mono_test = x_mono; 
            //     this->_nl_solver_master->max_it(50); 
            //     this->_nl_solver_master->solve(level.fun_monolithic(), x_mono_test); 

            //     // from monolithic to staggered ... 

            // }


            // level.transfer_from_stag1_to_monolithic(x_stag1, x_mono); 
            // level.transfer_from_stag2_to_monolithic(x_stag2, x_mono); 




            return true;
        }


        virtual void set_parameters(const Parameters params)
        {
            // AlternateMinimization<Matrix, Vector>::set_parameters(params);
            NonlinearSolver::set_parameters(params);
            TrustRegionBase::set_parameters(params);
        }






        bool solve(Function<Matrix, Vector> &fun, Vector &x_k) override 
        {
            std::cerr<<"-------AutoAdaptiveAlternateMin::solve not implemented ---------  \n"; 
            return false; 
        }





    private: 

        std::shared_ptr<TRSubproblem>        _stag1_subproblem;     
        std::shared_ptr<TRSubproblem>        _stag2_subproblem;   
        std::shared_ptr<TRSubproblem>        _monolithic_subproblem;   




    };

}
#endif //UTOPIA_AUTO_ADAPTIVE_ALTERNATE_MINIMIZATION_HPP
