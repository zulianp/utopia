/*
* @Author: alenakopanicakova
* @Date:   2016-08-10
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-12
*/
#ifndef UTOPIA_Linear_SOLVER_H
#define UTOPIA_Linear_SOLVER_H

#include <string>
#include "utopia_Core.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_Traits.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Monitor.hpp"


namespace  utopia 
{
    /**
     * @brief      The base class for linear solvers.
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class IterativeSolver : public LinearSolver<Matrix, Vector>,
                            public Monitor<Matrix, Vector>
    {
    
    public:
        typedef UTOPIA_SCALAR(Matrix)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Matrix)        SizeType;

        IterativeSolver(const Parameters params = Parameters())
        {
            set_parameters(params); 
        }
        
        virtual ~IterativeSolver( ){}
        
        /**
         * @brief      Solve routine. Needs to be provided by each solver.
         *
         * @param[in]  A     
         * @param[in]  b     
         * @param      x0    
         *
         * @return    
         */
        // virtual bool solve(const Matrix &A, const Vector &b, Vector &x0) = 0;

        inline void copy_parameters_from(const IterativeSolver &other) 
        {
            Parameters params;
            other.get_parameters(params);
            set_parameters(params);
        }

        virtual void get_parameters(Parameters &params) const
        {
            params.ksp_atol(atol_);
            params.ksp_rtol(rtol_);
            params.ksp_dtol(stol_);

            params.ksp_max_it(max_it_);
            params.linear_solver_verbose(verbose_);
            params.time_statistics(time_statistics_);
            params.linear_solver_time_statistics(time_statistics_);

            params.log_system(log_system_);
            params.log_iterates(log_iterates_);
        }


        virtual void set_parameters(const Parameters params) override
        {
            atol_               = params.ksp_atol();            
            rtol_               = params.ksp_rtol(); 
            stol_               = params.ksp_dtol(); 

            max_it_             = params.ksp_max_it(); 
            verbose_            = params.linear_solver_verbose(); 
            time_statistics_    = params.time_statistics();  
            time_statistics_    = params.linear_solver_time_statistics(); 

            log_system_         = params.log_system();
            log_iterates_       = params.log_iterates(); 
        }

        virtual bool apply(const Vector &rhs, Vector &sol) override
        {
            return this->solve(*this->get_operator(), rhs, sol);
        }
    
    protected:

        /**
         * @brief      Monitors(creating matlab script) iterate, hessian on given iterate.
         */
        virtual bool solver_monitor(const SizeType &it, Vector &x, Matrix &H) override
        {
            if(log_iterates_)
                monitor(it, x);
            
            if(log_system_)
                monitor(it, H); 

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
            {
                PrintInfo::print_init(method, status_variables); 
            }
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
            
            num_it_ = num_it; 
            conv_reason_ = convergence_reason; 

            if(verbose_ && mpi_world_rank() == 0)
            {
              ConvergenceReason::exitMessage(num_it, convergence_reason);
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
        virtual bool check_convergence(const SizeType &it, const Scalar &g_norm, const Scalar & r_norm, const Scalar &s_norm) override
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
            if( it >= max_it_)
            {
                exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
                return true; 
            }

            return false; 
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public:
        Scalar      atol() const               { return atol_; } 
        Scalar      rtol()  const              { return rtol_; } 
        Scalar      stol()  const              { return stol_; } 
        SizeType    max_it()  const            { return max_it_; } 

        bool      precondition() const          { return precondition_; } 
        bool      time_statistics() const       { return time_statistics_; } 


        bool log_iterates() const                { return log_iterates_; } 
        bool log_system() const                  { return log_system_; } 
        bool verbose() const                     { return verbose_; } 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar & atol_in ) { atol_ = atol_in; }; 
        void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; }; 
        void stol(const Scalar & stol_in ) { stol_ = stol_in; }; 
        void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; }; 
        void verbose(const bool & verbose_in ) {verbose_ = verbose_in; }; 
        
        void precondition(const bool & precondition_in ) { precondition_ = precondition_in; }; 
        void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; }; 
        void log_iterates(const bool & log_iterates_in ) { log_iterates_ = log_iterates_in; }; 
        void log_system(const bool & log_system_in ) { log_system_ = log_system_in; }; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        SizeType get_convergence_reason(){ return conv_reason_;  }
        SizeType get_num_it() { return num_it_;  }

    private:
        
        // ... GENERAL Iterative SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */  
        Scalar rtol_;                   /*!< Relative tolerance. */  
        Scalar stol_;                   /*!< Step tolerance. */  

        SizeType max_it_;               /*!< Maximum number of iterations. */  
        bool verbose_;                  /*!< Verobse enable? . */  

        bool time_statistics_;          /*!< Perform time stats or not? */  
        bool precondition_;             // todo


        bool log_iterates_;
        bool log_system_; 

        // to be passed out 
        SizeType num_it_; 
        SizeType  conv_reason_; 


        Chrono _time;                 /*!<Timing of solver. */
    };
}

#endif //UTOPIA_Linear_SOLVER_H
