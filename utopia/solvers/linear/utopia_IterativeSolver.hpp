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
                            public Monitor<Vector>
    {
    
    public:
        typedef UTOPIA_SCALAR(Matrix)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Matrix)        SizeType;

        IterativeSolver():  atol_(1e-9), rtol_(1e-9), stol_(1e-11), max_it_(300), verbose_(false), 
                            time_statistics_(true), log_iterates_(false), log_system_(false)
        {
            
        }
        
        virtual ~IterativeSolver( ){}
        
        virtual void read(Input &is) override
        {
            is.get("atol", atol_);
            is.get("rtol", rtol_);
            is.get("stol", stol_);

            is.get("max-it", max_it_);
            is.get("verbose", verbose_);
            is.get("time-statistics", time_statistics_);
            is.get("log-system", log_system_);
            is.get("log-iterates", log_iterates_);
        }

        virtual void print_usage(std::ostream &os) const override
        {
            os << "atol             : <real>\n";
            os << "rtol             : <real>\n";
            os << "stol             : <real>\n";
            os << "max-it           : <int>\n";
            os << "verbose          : <bool>\n";
            os << "time-statistics  : <bool>\n";
            os << "log-system       : <bool>\n";
            os << "log-iterates     : <bool>\n";
        }

        virtual bool apply(const Vector &rhs, Vector &sol) override
        {
            return this->solve(*this->get_operator(), rhs, sol);
        }
    
        Scalar get_time() { return _time.get_seconds();  }

    protected:

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
         * @brief      Writes CSV file with iteration info 
         *
         * @param[in]  it_global  The iterator global
         */
        virtual void print_statistics(const SizeType & it_global)
        {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if(!non_data_path.empty())
            {
                CSVWriter writer;
                if (mpi_world_rank() == 0)
                {
                    if(!writer.file_exists(non_data_path))
                    {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"num_its", "time"});
                    }
                    else
                        writer.open_file(non_data_path);
                    
                    writer.write_table_row<Scalar>({Scalar(it_global), this->get_time()});
                    writer.close_file();
                }
            }
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
        virtual Scalar      atol() const               { return atol_; } 
        virtual Scalar      rtol()  const              { return rtol_; } 
        virtual Scalar      stol()  const              { return stol_; } 

        virtual SizeType    max_it()  const            { return max_it_; } 

        virtual bool      time_statistics() const       { return time_statistics_; } 

        virtual bool log_iterates() const                { return log_iterates_; } 
        virtual bool log_system() const                  { return log_system_; } 
        virtual bool verbose() const                     { return verbose_; } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        virtual void atol(const Scalar & atol_in ) { atol_ = atol_in; }; 
        virtual void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; }; 
        virtual void stol(const Scalar & stol_in ) { stol_ = stol_in; }; 
        virtual void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; }; 
        virtual void verbose(const bool & verbose_in ) {verbose_ = verbose_in; }; 
        
        virtual void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; }; 
        virtual void log_iterates(const bool & log_iterates_in ) { log_iterates_ = log_iterates_in; }; 
        virtual void log_system(const bool & log_system_in ) { log_system_ = log_system_in; }; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        virtual SizeType get_convergence_reason(){ return conv_reason_;  }
        virtual SizeType get_num_it() { return num_it_;  }

    private:
        
        // ... GENERAL Iterative SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */  
        Scalar rtol_;                   /*!< Relative tolerance. */  
        Scalar stol_;                   /*!< Step tolerance. */  

        SizeType max_it_;               /*!< Maximum number of iterations. */  
        bool verbose_;                  /*!< Verobse enable? . */  

        bool time_statistics_;          /*!< Perform time stats or not? */  

        bool log_iterates_;
        bool log_system_; 

        // to be passed out 
        SizeType num_it_; 
        SizeType  conv_reason_; 


        Chrono _time;                 /*!<Timing of solver. */
    };
}

#endif //UTOPIA_Linear_SOLVER_H
