#ifndef UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP
#define UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP

#include "utopia_Core.hpp"


namespace utopia
{
    template<class Vector>
    class FunctionOperator final: public Operator<Vector> 
    {
        public:
            FunctionOperator(const std::function< void(const Vector &, Vector &) > operator_action)
            : operator_action_(operator_action)
            {}

            bool apply(const Vector &rhs, Vector &ret) const override
            {
                operator_action_(rhs, ret); 
                return true;
            }


        private:
            std::function< void(const Vector &, Vector &) > operator_action_; 
    };


    template<class Vector>
    class FunctionPreconditioner final: public Preconditioner<Vector> 
    {
        public:
            FunctionPreconditioner(const std::function< void(const Vector &, Vector &) > operator_action)
            : operator_action_(operator_action)
            {}

            bool apply(const Vector &rhs, Vector &ret) override
            {
                operator_action_(rhs, ret); 
                return true;
            }


        private:
            std::function< void(const Vector &, Vector &) > operator_action_; 
    };

    
    // TODO:: fix me, use precond with inverse action instead... 
    template<class Vector>
    class QuasiLinearSolver: public MatrixFreeLinearSolver<Vector>
    {
        public:
            virtual ~QuasiLinearSolver() {}

            virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol)
            {
                A.apply(rhs, sol); 
                return true; 
            }

    };



}

#endif //UTOPIA_MATRIX_FREE_SOLVER_INTERFACE_HPP