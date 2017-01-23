#ifndef UTOPIA_UTOPIA_MONITOR_HPP_HPP
#define UTOPIA_UTOPIA_MONITOR_HPP_HPP

namespace utopia {

    template<class Matrix, class Vector>
    class Monitor {
    public:
        typedef UTOPIA_SCALAR(Matrix)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Matrix) SizeType;

        virtual ~Monitor() {}

    public:

        /**
        * @brief      Monitors(creating matlab script) iterate, hessian on given iterate.
        */
        //FIXME move Matrix before Vector
        virtual bool solver_monitor(const SizeType &it, Vector &rhs, Matrix &system) = 0;

        /**
         * @brief      Initialization of nonlinear solver. Includes nice printout and starts calculating time of solve process.
         *
         * @param[in]  method            The method.
         * @param[in]  status_variables  The status variables.
         */
        virtual void init_solver(const std::string & method, const std::vector<std::string> status_variables) = 0;

        /**
         * @brief      Exit of solver.
         *
         * @param[in]  it              The number iterator
         * @param[in]  convergence_reason  The convergence reason
         */
        virtual void exit_solver(const SizeType &it, const Scalar & convergence_reason) = 0;

        /**
         * @brief      General function to check convergence in nonlinear solvers. It checks absolute, relative norm of gradient
         *             and lenght of the step size.
         *
         * @param[in]  norm_grad  The norm of the gradient.
         * @param[in]  rel_norm_grad  The relative norm of the gradient.
         * @param[in]  norm_step  The size of step.
         * @param[in]  it      The number of iterations.
         */
        virtual bool check_convergence(const SizeType &it, const Scalar & norm_grad, const Scalar &rel_norm_grad, const Scalar &norm_step) = 0;
    };
}

#endif //UTOPIA_UTOPIA_MONITOR_HPP_HPP
