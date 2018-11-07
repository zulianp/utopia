#ifndef UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
#define UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP

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

    template<class Matrix, class Vector>
    class NewtonBasedNonLinearSolver : public NonLinearSolver<Matrix, Vector>
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;


        NewtonBasedNonLinearSolver(const std::shared_ptr<Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(),
                        const Parameters &params = Parameters()): 
                        NonLinearSolver<Matrix, Vector>(params), 
                        linear_solver_(linear_solver)
        {
            set_parameters(params);        
        }

        virtual ~NewtonBasedNonLinearSolver() {}


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
         * @brief      Settter the parameters.
         *
         * @param[in]  params  The parameters
         */
        virtual void set_parameters(const Parameters params)
        {
            NonLinearSolver<Matrix, Vector>::set_parameters(params); 
            check_diff_         = params.differentiation_control(); 

        }


        /**
         * @brief      Changes linear solver used inside of nonlinear-solver. 
         *
         * @param[in]  linear_solver  The linear solver
         */
        virtual void set_linear_solver(const std::shared_ptr<Solver> &linear_solver)
        {
            linear_solver_ = linear_solver; 
        }

        inline DiffController &controller() { return controller_; }


    public:
        inline std::shared_ptr<Solver> linear_solver() const
        {
            return linear_solver_;
        }

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
        DiffController          controller_;
        bool                    check_diff_;        /*!< Enable differentiation control. */  

    };
}

#endif //UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
