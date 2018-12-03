#ifndef UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
#define UTOPIA_UTOPIA_NONLINEARSOLVER_HPP

#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"

#include "utopia_Parameters.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_ConjugateGradient.hpp"


namespace utopia 
{
    /**
     * @brief      The base class for all nonlinear solvers. Class provides basic functions used in all nonlinear solvers. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Vector>
    class NonLinearSolver : public Monitor<Vector>, public Configurable
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        NonLinearSolver(const Parameters &params = Parameters()): 
                        params_(params)
        {
            set_parameters(params);        
        }

        virtual ~NonLinearSolver() {}

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
        }


        virtual void read(Input &in) override
        {
            in.get("atol", atol_);
            in.get("rtol", rtol_);
            in.get("stol", stol_);
            in.get("max-it", max_it_);
            in.get("verbose", verbose_);
            in.get("time-statistics", time_statistics_);
        }

        virtual void print_usage(std::ostream &os) const override
        {
            os << "atol             : <real>\n";
            os << "rtol             : <real>\n";
            os << "stol             : <real>\n";
            os << "max-it           : <int>\n";
            os << "verbose          : <bool>\n";
            os << "time-statistics  : <bool>\n";
        }


protected:
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

            if(mpi_world_rank() == 0 && verbose_)
            {
                ConvergenceReason::exitMessage_nonlinear(num_it, convergence_reason);
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar & atol_in ) { atol_ = atol_in; }; 
        void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; }; 
        void stol(const Scalar & stol_in ) { stol_ = stol_in; }; 
        void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; }; 
        void verbose(const bool & verbose_in ) { verbose_ = verbose_in; }; 
        void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; }; 

        Scalar get_time() { return _time.get_seconds();  }


    protected:
        Parameters params_;                         /*!< Solver parameters. */  

        // ... GENERAL SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */  
        Scalar rtol_;                   /*!< Relative tolerance. */  
        Scalar stol_;                   /*!< Step tolerance. */  

        SizeType max_it_;               /*!< Maximum number of iterations. */  
        bool verbose_;              /*!< Verobse enable? . */  
        SizeType time_statistics_;      /*!< Perform time stats or not? */  

        Chrono _time;                 /*!<Timing of solver. */

    };



    template<class Vector>
    class MatrixFreeNonLinearSolver : public NonLinearSolver<Vector>
    {
    
    public:
        MatrixFreeNonLinearSolver(const Parameters &params = Parameters()):  
                                  NonLinearSolver<Vector>(params)
        {

        }

        virtual ~MatrixFreeNonLinearSolver() {}

        virtual bool solve(FunctionBase<Vector> &fun, Vector &x) = 0;

    };


}

#endif //UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
