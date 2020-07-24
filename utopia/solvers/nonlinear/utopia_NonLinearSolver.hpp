#ifndef UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
#define UTOPIA_UTOPIA_NONLINEARSOLVER_HPP

#include "utopia_ConvergenceReason.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_SolutionStatus.hpp"

namespace utopia {
    /**
     * @brief      The base class for all nonlinear solvers. Class provides basic functions used in all nonlinear
     * solvers.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Vector>
    class NonLinearSolver : public Monitor<Vector>, virtual public Configurable {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        NonLinearSolver() : atol_(1e-7), rtol_(1e-8), stol_(1e-9), max_it_(300), time_statistics_(true) {}

        ~NonLinearSolver() override = default;

        void read(Input &in) override {
            in.get("atol", atol_);
            in.get("rtol", rtol_);
            in.get("stol", stol_);

            in.get("max-it", max_it_);
            in.get("verbose", verbose_);
            in.get("time-statistics", time_statistics_);
        }

        void print_usage(std::ostream &os) const override {
            this->print_param_usage(os, "atol", "real", "Absolute tolerance.", std::to_string(atol()));
            this->print_param_usage(os, "rtol", "real", "Relative tolerance.", std::to_string(rtol()));
            this->print_param_usage(os, "stol", "real", "Step size tolerance.", std::to_string(stol()));

            this->print_param_usage(os, "max-it", "int", "Maximum number of iterations.", std::to_string(max_it()));
            this->print_param_usage(os, "verbose", "bool", "Turn on/off output.", std::to_string(verbose_));
            this->print_param_usage(
                os, "time-statistics", "bool", "Collect time-statistics.", std::to_string(time_statistics_));
        }

    protected:
        virtual void print_statistics(const SizeType &it_global) {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if (!non_data_path.empty()) {
                CSVWriter writer{};
                if (mpi_world_rank() == 0) {
                    if (!writer.file_exists(non_data_path)) {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"num_its", "time"});
                    } else
                        writer.open_file(non_data_path);

                    writer.write_table_row<Scalar>({Scalar(it_global), this->get_time()});
                    writer.close_file();
                }
            }
        }

        /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and starts calculating time of solve
         * process.
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        void init_solver(const std::string &method, const std::vector<std::string> status_variables) override {
            if (mpi_world_rank() == 0 && verbose_) {
                this->print_init_message(method, status_variables);
            }

            this->solution_status_.clear();
            _time.start();
        }

        virtual void print_init_message(const std::string &method,
                                        const std::vector<std::string> status_variables)  // override
        {
            if (mpi_world_rank() == 0 && verbose_) {
                PrintInfo::print_init(method, status_variables);
            }
        }

        /**
         * @brief      Exit of solver.
         *
         * @param[in]  num_it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        void exit_solver(const SizeType &num_it, const Scalar &convergence_reason) override {
            _time.stop();

            if (mpi_world_rank() == 0 && verbose_) {
                ConvergenceReason::exitMessage_nonlinear(num_it, convergence_reason);
                utopia::out() << "  Walltime of solve: " << _time.get_seconds() << " seconds. \n";
            }

            this->solution_status_.execution_time = _time.get_seconds();
            this->solution_status_.iterates = num_it;
            this->solution_status_.reason = convergence_reason;
        }

        /**
         * @brief      General function to check convergence in nonlinear solvers. It checks absolute, relative norm of
         * gradient and lenght of the step size.
         *
         * @param[in]  g_norm  The norm of the gradient.
         * @param[in]  r_norm  The relative norm of the gradient.
         * @param[in]  s_norm  The size of step.
         * @param[in]  it      The number of iterations.
         */
        bool check_convergence(const SizeType &it,
                               const Scalar &g_norm,
                               const Scalar &r_norm,
                               const Scalar &s_norm) override {
            bool converged = false;

            // termination because norm of grad is down
            if (g_norm < atol_) {
                exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
                this->solution_status_.reason = ConvergenceReason::CONVERGED_FNORM_ABS;
                converged = true;
            }

            // step size so small that we rather exit than wait for nan's
            else if (s_norm < stol_) {
                exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
                this->solution_status_.reason = ConvergenceReason::CONVERGED_SNORM_RELATIVE;
                converged = true;
            }

            // step size so small that we rather exit than wait for nan's
            else if (r_norm < rtol_) {
                exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
                this->solution_status_.reason = ConvergenceReason::CONVERGED_FNORM_RELATIVE;
                converged = true;
            }

            // check number of iterations
            else if (it > max_it_) {
                exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
                this->solution_status_.reason = ConvergenceReason::DIVERGED_MAX_IT;
                converged = true;
            }

            if (converged) {
                this->solution_status_.iterates = it;
                this->solution_status_.gradient_norm = g_norm;
                this->solution_status_.relative_gradient_norm = r_norm;
                this->solution_status_.step_norm = s_norm;
            }

            return converged;
        }

    public:
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Scalar atol() const { return atol_; }
        Scalar rtol() const { return rtol_; }
        Scalar stol() const { return stol_; }
        SizeType max_it() const { return max_it_; }
        bool verbose() const { return verbose_; }
        bool time_statistics() const { return time_statistics_; }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void atol(const Scalar &atol_in) { atol_ = atol_in; };
        void rtol(const Scalar &rtol_in) { rtol_ = rtol_in; };
        void stol(const Scalar &stol_in) { stol_ = stol_in; };
        void max_it(const SizeType &max_it_in) { max_it_ = max_it_in; };
        void verbose(const bool &verbose_in) { verbose_ = verbose_in; };
        void time_statistics(const bool &time_statistics_in) { time_statistics_ = time_statistics_in; };

        Scalar get_time() { return _time.get_seconds(); }

    protected:
        Scalar atol_; /*!< Absolute tolerance. */
        Scalar rtol_; /*!< Relative tolerance. */
        Scalar stol_; /*!< Step tolerance. */

        SizeType max_it_;          /*!< Maximum number of iterations. */
        bool verbose_{false};      /*!< Verobse enable? . */
        SizeType time_statistics_; /*!< Perform time stats or not? */

        Chrono _time; /*!<Timing of solver. */
    };

    template <class Vector>
    class MatrixFreeNonLinearSolver : public NonLinearSolver<Vector> {
    public:
        MatrixFreeNonLinearSolver() : NonLinearSolver<Vector>() {}

        ~MatrixFreeNonLinearSolver() override = default;

        virtual bool solve(FunctionBase<Vector> &fun, Vector &x) = 0;
    };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_NONLINEARSOLVER_HPP
