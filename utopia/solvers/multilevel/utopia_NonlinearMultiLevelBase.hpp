#ifndef UTOPIA_NONLINEAR_ML_BASE_HPP
#define UTOPIA_NONLINEAR_ML_BASE_HPP
#include "utopia_Level.hpp"
#include "utopia_Transfer.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_MatrixTransfer.hpp"


#include "utopia_MultiLevelEvaluations.hpp"



#include <algorithm>
#include <vector>

namespace utopia {
#define CHECK_NUM_PRECISION_mode


    /**
     * @brief      Base class for all nonlinear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy - calls into assembly routines on each level. \n
     *             Different levels are created by interpolation and restriction.\n
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class NonlinearMultiLevelBase : public MultiLevelBase<Matrix, Vector>, public Monitor<Matrix, Vector> {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Transfer<Matrix, Vector> Transfer;
        typedef utopia::ExtendedFunction<Matrix, Vector> Fun;
        typedef std::shared_ptr<Fun> FunPtr;
        typedef utopia::MatrixTransfer<Matrix, Vector> MatrixTransfer;

        using MultiLevelBase<Matrix, Vector>::set_transfer_operators;

        NonlinearMultiLevelBase(const Parameters params = Parameters())
        {
            set_parameters(params);
        }

        virtual ~NonlinearMultiLevelBase(){}


        /**
         * @brief      Sets the parameters.
         *
         * @param[in]  params  The parameters
         */
        virtual void set_parameters(const Parameters params) override
        {
            MultiLevelBase<Matrix, Vector>::set_parameters(params);

            atol_               = params.atol();
            rtol_               = params.rtol();
            stol_               = params.stol();

            max_it_             = params.max_it();
            verbose_            = params.verbose();
            time_statistics_    = params.time_statistics();
        }



        /**
         * @brief      The solve function for nonlinear multilevel solvers.
         *
         * @param[in]  rhs   Function to be minimized.
         * @param      x_h   The initial guess.
         *
         */
        virtual bool solve( Vector &x_h) = 0;


        /**
         * @brief      Fnction inits functions associated with assemble on each level.
         *
         * @param[in]  level_functions  The level functions
         *
         */
        virtual bool set_functions(const std::vector<FunPtr> &level_functions)
        {
            level_functions_.clear();
            level_functions_.insert(level_functions_.begin(), level_functions.begin(), level_functions.end());
            return true;
        }

        /* @brief
         Function initializes transfer  operators.
         Operators need to be ordered FROM COARSE TO FINE.
         *
         * @param[in]  interpolation_operators                The interpolation operators.
         * @param[in]  projection_operators                   The projection operators.
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators)
        {
            this->transfers_.clear();
            for(auto I = interpolation_operators.begin(), P = projection_operators.begin(); I != interpolation_operators.end() && P != projection_operators.end(); ++I, ++P )
                this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I, *P));

            return true;
        }

        /**
         * @brief      Sets the transfer operators.
         *
         * @param[in]  interpolation_operators  The interpolation operators
         * @param[in]  restriction_operators    The restriction operators
         * @param[in]  projection_operators     The projection operators
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &restriction_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators)
        {
            this->transfers_.clear();
            for(auto I = interpolation_operators.begin(), R = restriction_operators.begin(), P = projection_operators.begin(); I != interpolation_operators.end() && R != restriction_operators.end() &&  P != projection_operators.end(); ++I, ++R, ++P )
                this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I, *R, *P));

            return true;
        }


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Scalar      atol() const               { return atol_; }
        Scalar      rtol()  const              { return rtol_; }
        Scalar      stol()  const              { return stol_; }
        SizeType    max_it()  const            { return max_it_; }
        bool        verbose() const            { return verbose_; }
        bool        time_statistics() const    { return time_statistics_; }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar & atol_in ) { atol_ = atol_in; };
        void rtol(const Scalar & rtol_in ) { rtol_ = rtol_in; };
        void stol(const Scalar & stol_in ) { stol_ = stol_in; };
        void max_it(const SizeType & max_it_in ) { max_it_ = max_it_in; };
        void verbose(const bool & verbose_in ) { verbose_ = verbose_in; };
        void time_statistics(const bool & time_statistics_in ) { time_statistics_ = time_statistics_in; };


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
                PrintInfo::print_init(method, status_variables);

            _time.start();
        }


        virtual void print_init_message(const std::string &method, const std::vector<std::string> status_variables)
        {
            if(mpi_world_rank() == 0 && verbose_)
                PrintInfo::print_init(method, status_variables);
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

            status_.reason = convergence_reason;
            status_.iterates = num_it;

            if(verbose_)
            {
                ConvergenceReason::exitMessage_nonlinear(num_it, convergence_reason);
                if(mpi_world_rank() == 0)
                    std::cout<<"  Walltime of solve: " << _time.get_seconds() << " seconds. \n";
            }
        }


        virtual bool solver_monitor(const SizeType& /*it*/, Vector & /*x*/, Matrix & /*H*/) override
        {
            std::cout<<"utopia::NonlinearMultilevelBase:: WE ARE NOT SUPPORTING this function at the moment... \n";
            return true;
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


        /**
         * @brief      Function looks up for ids, where we should apply Dirichlet BC and set value to required one
         *
         * @param      fun   The fun
         * @param      x
         *
         */
        virtual bool make_iterate_feasible(Fun & fun, Vector & x)
        {
            Vector bc_values;
            fun.get_eq_constrains_values(bc_values);

            Vector bc_ids;
            fun.get_eq_constrains_flg(bc_ids);

            if(local_size(bc_ids).get(0) != local_size(bc_values).get(0)) {
                std::cerr<<"utopia::NonlinearMultiLevelBase::make_iterate_feasible:: local sizes do not match... \n";
            }

            {
                Write<Vector> w(x);
                Read<Vector>  r_id(bc_ids);
                Read<Vector>  r_val(bc_values);

                Range range_w = range(x);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    Scalar id = bc_ids.get(i);
                    Scalar value = bc_values.get(i);

                    if(id == 1)
                        x.set(i, value);
                }
            }

            return true;
        }


        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      c     The correction
         */
        virtual bool zero_correction_related_to_equality_constrain(Fun & fun, Vector & c)
        {
            Vector bc;
            fun.get_eq_constrains_flg(bc);

            {
                Write<Vector> w(c);
                Read<Vector> r(bc);

                Range range_w = range(c);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(bc.get(i) == 1)
                        c.set(i, 0);
                }
            }
            return true;
        }


        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      M     matrix
         *
         */
        virtual bool zero_correction_related_to_equality_constrain_mat(Fun & fun, Matrix & M)
        {
            Vector bc;
            fun.get_eq_constrains_flg(bc);
            std::vector<SizeType> index;

            {
                Read<Vector> r(bc);

                Range range_w = range(bc);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(bc.get(i) == 1)
                        index.push_back(i);
                }
            }
            set_zero_rows(M, index);

            // horible solution....
            {
                ReadAndWrite<Matrix> w(M);
                Range r = row_range(M);

                //You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
                for(SizeType i = r.begin(); i != r.end(); ++i)
                {
                    if(std::abs(M.get(i,i)) < 1e-15)
                        M.set(i, i, 1.0);
                }
            }


            return true;
        }



        /**
         * @return     Name of solver - to have nice printouts
         */
        virtual std::string name() = 0;

        /**
         * @brief      Init internal memory used for implementation of given multilevel solver
         *
         * @param[in]  fine_local_size  The local size of fine level problem
         */
        virtual void init_memory(const SizeType & fine_local_size) = 0;



        inline Fun &function(const SizeType level)
        {
            return *level_functions_[level];
        }

        inline const Fun &function(const SizeType level) const
        {
            return *level_functions_[level];
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
                        writer.write_table_row<std::string>({"v_cycles", "time"});
                    }
                    else
                        writer.open_file(non_data_path);

                    writer.write_table_row<Scalar>({Scalar(it_global), this->get_time()});
                    writer.close_file();
                }
            }
        }



    protected:
        std::vector<FunPtr>                      level_functions_;

        // ... GENERAL SOLVER PARAMETERS ...
        Scalar atol_;                   /*!< Absolute tolerance. */
        Scalar rtol_;                   /*!< Relative tolerance. */
        Scalar stol_;                   /*!< Step tolerance. */

        SizeType max_it_;               /*!< Maximum number of iterations. */
        SizeType verbose_;              /*!< Verobse enable? . */
        SizeType time_statistics_;      /*!< Perform time stats or not? */


        Chrono _time;                 /*!<Timing of solver. */

        SolutionStatus status_;

    };

}

#endif //UTOPIA_NONLINEAR_ML_BASE_HPP

