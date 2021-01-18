#ifndef UTOPIA_UTOPIA_MONITOR_HPP_HPP
#define UTOPIA_UTOPIA_MONITOR_HPP_HPP

#include "utopia_SolutionStatus.hpp"

namespace utopia {

    template <class Vector>
    class Monitor {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        virtual ~Monitor() = default;

    public:
        /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and
         * starts calculating time of solve process.
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        virtual void init_solver(const std::string &method, const std::vector<std::string> status_variables) = 0;

        /**
         * @brief      Exit of solver.
         *
         * @param[in]  it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        virtual void exit_solver(const SizeType &it, const Scalar &convergence_reason) = 0;

        /**
         * @brief      General function to check convergence in nonlinear solvers. It
         * checks absolute, relative norm of gradient and lenght of the step size.
         *
         * @param[in]  norm_grad  The norm of the gradient.
         * @param[in]  rel_norm_grad  The relative norm of the gradient.
         * @param[in]  norm_step  The size of step.
         * @param[in]  it      The number of iterations.
         */
        virtual bool check_convergence(const SizeType &it,
                                       const Scalar &norm_grad,
                                       const Scalar &rel_norm_grad,
                                       const Scalar &norm_step) = 0;

        const SolutionStatus &solution_status() const { return solution_status_; }

        void solution_status(const SolutionStatus &sol) { solution_status_ = sol; }

        SizeType get_num_it() const { return solution_status_.iterates; }

        SizeType get_convergence_reason() const { return solution_status_.reason; }

    protected:
        SolutionStatus solution_status_;
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_MONITOR_HPP_HPP
