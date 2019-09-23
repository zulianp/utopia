#ifndef UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
#define UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP

#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_NonLinearSolver.hpp"

namespace utopia {


    enum InexactNewtonForcingStartegies  {  ZERO            = 0,  
                                            CAI             = 1,
                                            DEMBO           = 2,
                                            SUPERLINEAR     = 3,
                                            QUADRATIC       = 4,
                                            QUADRATIC_2     = 5 };   


    template<class Matrix, class Vector>
    class NewtonBase : public NonLinearSolver<Vector> {
    public:
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;
        typedef utopia::LinearSolver<Matrix, Vector>    Solver;

        NewtonBase(const std::shared_ptr<Solver> &linear_solver)
        : NonLinearSolver<Vector>(), linear_solver_(linear_solver), check_diff_(false), forcing_strategy_(InexactNewtonForcingStartegies::ZERO)
        { }

        virtual ~NewtonBase() {}


        virtual bool solve(Function<Matrix, Vector> &fun, Vector &x) = 0;

        virtual bool solve(ExtendedFunction<Matrix, Vector> &fun, Vector &x, const Vector & rhs)
        {
            fun.set_rhs(rhs);
            bool converged = this->solve(fun, x);
            fun.reset_rhs();
            return converged;
        }

        virtual void forcing_strategy(const InexactNewtonForcingStartegies & strategy)
        {
            forcing_strategy_ = strategy; 
        }

        virtual InexactNewtonForcingStartegies forcing_strategy()
        {
            return forcing_strategy_; 
        }    

        virtual bool has_forcing_strategy()
        {
            return (forcing_strategy_>0) ? true : false; 
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

        virtual void read(Input &in) override
        {
            NonLinearSolver<Vector>::read(in);
            in.get("check_diff", check_diff_);

            if(linear_solver_) {
                in.get("linear-solver", *linear_solver_);
            }
        }

        virtual void print_usage(std::ostream &os) const override
        {
            NonLinearSolver<Vector>::print_usage(os);
            this->print_param_usage(os, "check_diff", "bool", "Enables finite difference controller", std::to_string(check_diff_));

            if(linear_solver_) {
                this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "-");
                linear_solver_->print_usage(os);
            } else {
                this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for linear solver.", "- (null)");
            }
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
            this->solution_status_.num_linear_solves++;
            return linear_solver_->apply(rhs, sol);
        }

        inline bool has_preconditioned_solver()
        {
            return dynamic_cast< PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get());
        }


        inline bool linear_solve(const Matrix &mat, const Matrix &prec, const Vector &rhs, Vector &sol)
        {
            static_cast< PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get())->update(make_ref(mat), make_ref(prec));
            this->solution_status_.num_linear_solves++;
            return linear_solver_->apply(rhs, sol);
        }


        virtual Scalar estimate_ls_atol(const Scalar & gnorm, const Scalar & it)
        {
            if(forcing_strategy_==InexactNewtonForcingStartegies::CAI)
            {
                return 10e-4*gnorm; 
            }   
            else if(forcing_strategy_ == InexactNewtonForcingStartegies::DEMBO)
            {
                return (gnorm * std::min(std::sqrt(gnorm), 1./(it+1))); 
            }
            else if(forcing_strategy_ == InexactNewtonForcingStartegies::SUPERLINEAR)
            {
                return (gnorm *  std::min(std::sqrt(gnorm), 0.5)); 
            }
            else if(forcing_strategy_ == InexactNewtonForcingStartegies::QUADRATIC)
            {
                return (gnorm * std::min(gnorm, 1./(it+1))); 
            }                        
            else if(forcing_strategy_ == InexactNewtonForcingStartegies::QUADRATIC_2)
            {
                return (gnorm * std::min(gnorm, 0.5)); 
            }               
            else
            {
                utopia_error("utopia::Newton:: invalid choice of forcing strategy.... \n"); 
                return 0.0; 
            }

            return 0.0; 

        }   



        std::shared_ptr<Solver> linear_solver_;     /*!< Linear solver parameters. */
        DiffController          controller_;
        bool                    check_diff_;        /*!< Enable differentiation control. */

        InexactNewtonForcingStartegies forcing_strategy_; 

    };
}

#endif //UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
