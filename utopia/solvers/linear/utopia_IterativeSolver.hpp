#ifndef UTOPIA_Linear_SOLVER_H
#define UTOPIA_Linear_SOLVER_H

#include <string>
#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    /**
     * @brief      The base class for linear solvers.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class IterativeSolver : public LinearSolver<Matrix, Vector>, virtual public Monitor<Vector> {
    public:
        using Scalar = typename utopia::Traits<Matrix>::Scalar;
        using SizeType = typename utopia::Traits<Matrix>::SizeType;

        IterativeSolver() = default;

        // IterativeSolver(const IterativeSolver &other) = default;
        IterativeSolver(IterativeSolver &&other) = default;
        IterativeSolver &operator=(const IterativeSolver &) = default;
        IterativeSolver &operator=(IterativeSolver &&) = default;

        ~IterativeSolver() override = default;

        void read(Input &in) override {
            LinearSolver<Matrix, Vector>::read(in);

            in.get("atol", atol_);
            in.get("rtol", rtol_);
            in.get("stol", stol_);

            in.get("max_it", max_it_);
            in.get_deprecated("max-it", "max_it", max_it_);
        }

        void print_usage(std::ostream &os) const override {
            LinearSolver<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "atol", "real", "Absolute tolerance.", std::to_string(atol_init_));
            this->print_param_usage(os, "rtol", "real", "Relative tolerance.", std::to_string(rtol_init_));
            this->print_param_usage(os, "stol", "real", "Minimum step-size.", std::to_string(stol_init_));
            this->print_param_usage(os, "max_it", "int", "Maximum number of iterations.", std::to_string(max_it_init_));
        }

        bool apply(const Vector &rhs, Vector &sol) override { return this->solve(*this->get_operator(), rhs, sol); }

        virtual bool smooth(const Vector &rhs, Vector &x) {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            this->solve(*this->get_operator(), rhs, x);
            this->max_it(temp);
            return true;
        }

        SizeType sweeps() { return max_it(); }

        void sweeps(const SizeType &sweeps) { return max_it(sweeps); }

        Scalar get_time() { return _time.get_seconds(); }

        IterativeSolver<Matrix, Vector> *clone() const override = 0;

        IterativeSolver<Matrix, Vector>(const IterativeSolver<Matrix, Vector> &other)
            : LinearSolver<Matrix, Vector>(other),
              atol_(other.atol_),
              rtol_(other.rtol_),
              stol_(other.stol_),
              max_it_(other.max_it_),
              norm_freq_(other.norm_freq_) {}

    protected:
        /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and
         * starts calculating time of solve process.
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        void init_solver(const std::string &method, const std::vector<std::string> status_variables) override {
            if (mpi_world_rank() == 0 && this->verbose()) {
                PrintInfo::print_init(method, status_variables);
            }

            this->solution_status_.clear();
            _time.start();
        }

        /**
         * @brief      Writes CSV file with iteration info
         *
         * @param[in]  it_global  The iterator global
         */
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
         * @brief      Exit of solver.
         *
         * @param[in]  num_it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        void exit_solver(const SizeType &num_it, const Scalar &convergence_reason) override {
            _time.stop();

            if (mpi_world_rank() == 0 && this->verbose()) {
                ConvergenceReason::exitMessage(num_it, convergence_reason);
                if (mpi_world_rank() == 0)
                    utopia::out() << "  Walltime of solve: " << _time.get_seconds() << " seconds. \n";
            }

            this->solution_status_.execution_time = _time.get_seconds();
        }

        /**
         * @brief      General function to check convergence in nonlinear solvers. It
         * checks absolute, relative norm of gradient and lenght of the step size.
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
            if (compute_norm(it)) {
                // termination because norm of grad is down
                if (g_norm <= atol_) {
                    exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
                    this->solution_status_.reason = ConvergenceReason::CONVERGED_FNORM_ABS;
                    converged = true;
                }

                // step size so small that we rather exit than wait for nan's
                if (s_norm <= stol_) {
                    exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
                    this->solution_status_.reason = ConvergenceReason::CONVERGED_SNORM_RELATIVE;
                    converged = true;
                }

                // step size so small that we rather exit than wait for nan's
                if (r_norm <= rtol_) {
                    exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
                    this->solution_status_.reason = ConvergenceReason::CONVERGED_FNORM_RELATIVE;
                    converged = true;
                }
            }

            // check number of iterations
            if (it >= max_it_) {
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

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public:
        virtual Scalar atol() const { return atol_; }
        virtual Scalar rtol() const { return rtol_; }
        virtual Scalar stol() const { return stol_; }

        virtual SizeType max_it() const { return max_it_; }

        virtual SizeType norm_frequency() const { return norm_freq_; }

        virtual bool compute_norm(const SizeType &it) const {
            if (norm_freq_ == 0.0) {
                return false;
            } else {
                return (it % norm_freq_ == 0) ? true : false;
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        virtual void atol(const Scalar &atol_in) { atol_ = atol_in; };
        virtual void rtol(const Scalar &rtol_in) { rtol_ = rtol_in; };
        virtual void stol(const Scalar &stol_in) { stol_ = stol_in; };
        virtual void max_it(const SizeType &max_it_in) { max_it_ = max_it_in; };

        /**
         * @brief Define frequency with which we compute norms
         * @details 0 - never, any other number uses modulo
         *
         * @param freq - defines how often we compute norms
         */
        virtual void norm_frequency(const SizeType &freq) { norm_freq_ = freq; };

    private:
        static constexpr Scalar atol_init_{1e-9};
        static constexpr Scalar rtol_init_{1e-9};
        static constexpr Scalar stol_init_{1e-11};
        static constexpr SizeType max_it_init_{1000};

        // FIXME these fields should be removed and set directly in the backend state
        // variables
        // ... GENERAL Iterative SOLVER PARAMETERS ...
        Scalar atol_{atol_init_}; /*!< Absolute tolerance. */
        Scalar rtol_{rtol_init_}; /*!< Relative tolerance. */
        Scalar stol_{stol_init_}; /*!< Step tolerance. */

        SizeType max_it_{max_it_init_}; /*!< Maximum number of iterations. */

        Chrono _time; /*!<Timing of solver. */
        SizeType norm_freq_{1};
    };
}  // namespace utopia

#endif  // UTOPIA_Linear_SOLVER_H
