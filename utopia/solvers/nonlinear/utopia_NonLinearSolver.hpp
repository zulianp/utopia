#ifndef UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
#define UTOPIA_UTOPIA_NONLINEARSOLVER_HPP

#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"

#include "utopia_Parameters.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_PreconditionedSolver.hpp"


namespace utopia 
{
    /**
     * @brief      The base class for all nonlinear solvers. Class provides basic functions used in all nonlinear solvers. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class NonLinearSolver : public Monitor<Matrix, Vector>
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;


        NonLinearSolver(const std::shared_ptr<Solver> &linear_solver = std::shared_ptr<Solver>(),
                        const Parameters &params = Parameters())
        : linear_solver_(linear_solver),
          params_(params)
        {
            set_parameters(params);        
        }


        virtual ~NonLinearSolver() {}

        virtual bool solve(Function<Matrix, Vector> &fun, Vector &x) = 0;


        virtual bool solve(ExtendedFunction<Matrix, Vector> &fun, Vector &x, const Vector & rhs)
        {
            fun.set_rhs(rhs); 
            bool converged = this->solve(fun, x); 
            fun.reset_rhs(); 
            return converged; 
        }



        /**
         * @brief      Enables the differentiation control.
         *
         * @param[in]  checkDiff  Option, if eanable diff_control or no. 
         */
        void enable_differentiation_control(bool checkDiff) 
        {
            check_diff_ = checkDiff; 
        }

        inline bool differentiation_control_enabled() const 
        {
            return check_diff_; 
        }

        bool check_values(const SizeType iterations, const Function<Matrix, Vector> &fun, const Vector &x, const Vector &gradient, const Matrix &hessian)
        {
            if (check_diff_ && !controller_.check(fun, x, gradient, hessian)) 
            {
                exit_solver(iterations, norm2(gradient), ConvergenceReason::DIVERGED_INNER); 
                return false;
            }

            return true;
        }

        /**
         * @brief      Getter for parameters. 
         */
        Parameters parameters()
        {
            return params_;
        }

        /**
         * @brief      Settter the parameters.
         *
         * @param[in]  params  The parameters
         */
        virtual void set_parameters(const Parameters params)
        {
            atol_               = params.atol();            
            rtol_               = params.rtol(); 
            stol_               = params.stol(); 

            max_it_             = params.max_it(); 
            verbose_            = params.verbose(); 
            time_statistics_    = params.time_statistics();  

            log_iterates_       = params.log_iterates(); 
            log_system_         = params.log_system(); 
            check_diff_         = params.differentiation_control(); 

            linear_solver_->set_parameters(params); 
        }


        /**
         * @brief      Changes linear solver used inside of nonlinear-solver. 
         *
         * @param[in]  linear_solver  The linear solver
         */
        void set_linear_solver(const std::shared_ptr<Solver> &linear_solver = std::shared_ptr<Solver>())
        {
            linear_solver_ = linear_solver; 
        }



        inline DiffController &controller() { return controller_; }


protected:
        /**
         * @brief      Monitors(creating matlab script) iterate, hessian on given iterate.
         */
        virtual bool solver_monitor(const SizeType& it, Vector & x, Matrix & H) override
        {
            if(log_iterates_)
            {
                monitor(it, x);
            }
            if(log_system_)
            {
                monitor(it, H); 
            }

            return true; 
        }


        /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and starts calculating time of solve process. 
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        virtual void init_solver(const std::string &method, const std::vector<std::string> status_variables) override
        {
            if(mpi_world_rank() == 0 && verbose_)
                PrintInfo::print_init(method, status_variables); 
            
            _time.start();
        }     


        /**
         * @brief      Exit of solver. 
         *
         * @param[in]  num_it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        virtual void exit_solver(const SizeType &num_it, const Scalar & convergence_reason) override
         {            
            _time.stop();

            params_.convergence_reason(convergence_reason);
            params_.num_it(num_it);

            if(verbose_)
            {
                ConvergenceReason::exitMessage_nonlinear(num_it, convergence_reason);

                if(mpi_world_rank() == 0)
                    std::cout<<"  Walltime of solve: " << _time.get_seconds() << " seconds. \n";
                    
            }
         }



         /**
          * @brief      General function to check convergence in nonlinear solvers. It checks absolute, relative norm of gradient
          *             and lenght of the step size.   
          *
          * @param[in]  g_norm  The norm of the gradient. 
          * @param[in]  r_norm  The relative norm of the gradient. 
          * @param[in]  s_norm  The size of step. 
          * @param[in]  it      The number of iterations. 
          */
        virtual bool check_convergence(const SizeType &it, const Scalar & g_norm, const Scalar & r_norm, const Scalar & s_norm) override
        {   
            // termination because norm of grad is down
            if(g_norm < atol_)
            {
                exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
                return true; 
            }

            // step size so small that we rather exit than wait for nan's
            if(s_norm < stol_)
            {
                exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
                return true; 
            }

            // step size so small that we rather exit than wait for nan's
            if(r_norm < rtol_)
            {
                exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
                return true; 
            }

            // check number of iterations
            if( it > max_it_)
            {
                exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
                return true; 
            }

            return false; 
        }


public:
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Scalar      atol() const               { return atol_; } 
        Scalar      rtol()  const              { return rtol_; } 
        Scalar      stol()  const              { return stol_; } 
        SizeType    max_it()  const            { return max_it_; } 
        bool        verbose() const                     { return verbose_; } 
        bool        time_statistics() const       { return time_statistics_; } 
        bool        log_iterates() const          {return log_iterates_; }
        bool        log_system() const          {return log_system_; }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar & atol_in ) { atol_ = atol_in; }; 
        void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; }; 
        void stol(const Scalar & stol_in ) { stol_ = stol_in; }; 
        void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; }; 
        void verbose(const bool & verbose_in ) { verbose_ = verbose_in; }; 
        void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; }; 
        void log_iterates(const bool & log_iterates) { log_iterates_  = log_iterates; }; 
        void log_system(const bool & log_system) { log_system_  = log_system; }; 


        Scalar get_time() { return _time.get_seconds();  }

    protected:
        inline bool linear_solve(const Matrix &mat, const Vector &rhs, Vector &sol)
        {
            linear_solver_->update(make_ref(mat));
            return linear_solver_->apply(rhs, sol);
        }

        inline bool has_preconditioned_solver()
        {
            return dynamic_cast< PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get());
        }


        inline bool linear_solve(const Matrix &mat, const Matrix &prec, const Vector &rhs, Vector &sol)
        {
            static_cast< PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get())->update(make_ref(mat), make_ref(prec));
            return linear_solver_->apply(rhs, sol);
        }


        std::shared_ptr<Solver> linear_solver_;     /*!< Linear solver parameters. */  
        Parameters params_;        /*!< Solver parameters. */  
        DiffController controller_;

        // ... GENERAL SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */  
        Scalar rtol_;                   /*!< Relative tolerance. */  
        Scalar stol_;                   /*!< Step tolerance. */  

        SizeType max_it_;               /*!< Maximum number of iterations. */  
        SizeType verbose_;              /*!< Verobse enable? . */  
        SizeType time_statistics_;      /*!< Perform time stats or not? */  

        bool log_iterates_;             /*!< Monitoring of iterate. */  
        bool log_system_;               /*!< Monitoring of hessian/jacobian. */  
        bool check_diff_;               /*!< Enable differentiation control. */  


        Chrono _time;                 /*!<Timing of solver. */

    };
}

#endif //UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
