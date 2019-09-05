#ifndef UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
#define UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP


#include "utopia_Preconditioner.hpp"

namespace  utopia
{

    template<class Vector>
    class MatrixFreeLinearSolver : virtual public Configurable
    {
        public:
            virtual ~MatrixFreeLinearSolver() {}
            virtual bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) = 0;

            virtual MatrixFreeLinearSolver * clone() const =0;

            virtual void read(Input &/*in*/) override{ }
            virtual void print_usage(std::ostream & /*os*/) const override{ }
    };


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
    class EmptyPrecondMatrixFreeLinearSolver final: public MatrixFreeLinearSolver<Vector>
    {
        public:
            void set_preconditioner(const std::shared_ptr<Preconditioner<Vector> > &precond)
            {
                precond_ = precond;
            }

            bool solve(const Operator<Vector> &/*A*/, const Vector &rhs, Vector &sol) override
            {
                if(precond_){
                    precond_->apply(rhs, sol);
                }
                else{
                    utopia_warning("EmptyPrecondMatrixFreeLinearSolver: preconditioner is missing \n");
                }
                return true;
            }


            EmptyPrecondMatrixFreeLinearSolver * clone() const override
            {
                return new EmptyPrecondMatrixFreeLinearSolver(*this);
            }

            void read(Input &in) override
            {
                MatrixFreeLinearSolver<Vector>::read(in);
                if(precond_){
                    in.get("precond", *precond_);
                }
            }


            void print_usage(std::ostream &os) const override
            {
                MatrixFreeLinearSolver<Vector>::print_usage(os);
                this->print_param_usage(os, "precond", "Preconditioner", "Input parameters for preconditioner.", "-");
            }

        private:
            std::shared_ptr<Preconditioner<Vector> > precond_;

    };

}

#endif //UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
