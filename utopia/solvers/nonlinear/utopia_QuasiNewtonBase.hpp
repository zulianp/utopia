#ifndef UTOPIA_QUASI_NEWTON_BASE_HPP
#define UTOPIA_QUASI_NEWTON_BASE_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_HessianApproximations.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Vector>
    class QuasiNewtonBase
    {
        typedef UTOPIA_SCALAR(Vector)                               Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                            SizeType;
        
        typedef utopia::HessianApproximation<Vector>                HessianApproximation;
        typedef utopia::MatrixFreeLinearSolver<Vector>              MFSolver;
        
    public:

        QuasiNewtonBase(    const std::shared_ptr <HessianApproximation> &hessian_approx, 
                            const std::shared_ptr <MFSolver> &solver): 
                            hessian_approx_strategy_(hessian_approx), 
                            mf_linear_solver_(solver)
        {

        }
        

        virtual ~QuasiNewtonBase() {}
        

        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
        {
            hessian_approx_strategy_      = strategy;
            return true;
        }

        virtual void set_linear_solver(const std::shared_ptr<MFSolver> &linear_solver)
        {
            mf_linear_solver_ = linear_solver; 
        }
        
        inline virtual std::shared_ptr<MFSolver> linear_solver() const
        {
            return mf_linear_solver_;
        }

        virtual void update(const Vector & s, const Vector & y)
        {
            hessian_approx_strategy_->update(s, y);
        }

    protected:         
        inline bool linear_solve(const Vector &rhs, Vector &sol)
        {
            auto multiplication_action = FunctionOperator<Vector>(hessian_approx_strategy_->get_apply_H()); 
            return mf_linear_solver_->solve(multiplication_action, rhs, sol);             
        }

    protected:
        std::shared_ptr<HessianApproximation>   hessian_approx_strategy_;    
        std::shared_ptr<MFSolver>               mf_linear_solver_;   

    };
    
}
#endif //UTOPIA_QUASI_NEWTON_BASE_HPP
