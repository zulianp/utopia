#ifndef UTOPIA_QP_SOLVER_HPP
#define UTOPIA_QP_SOLVER_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"

#include <cmath>

namespace utopia
{

    template<class Matrix, class Vector>
    class QPSolver :    public virtual PreconditionedSolver<Matrix, Vector>,
                        public virtual VariableBoundSolverInterface<Vector> {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            QPSolver() {}
            virtual ~QPSolver() {}
            virtual QPSolver * clone() const override = 0;

            virtual void init_memory(const SizeType & ls)  override 
            {
                VariableBoundSolverInterface<Vector>::init_memory(ls); 
                PreconditionedSolver<Matrix, Vector>::init_memory(ls); 
            }      
    };


    template<class Vector>
    class MatrixFreeQPSolver : public virtual MatrixFreeLinearSolver<Vector>,
                               public virtual VariableBoundSolverInterface<Vector> {
        public:
            typedef UTOPIA_SCALAR(Vector)                           Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

            MatrixFreeQPSolver()
            {}

            virtual ~MatrixFreeQPSolver()
            {}

            virtual MatrixFreeQPSolver * clone() const override = 0;

            virtual void init_memory(const SizeType & ls)  override 
            {
                VariableBoundSolverInterface<Vector>::init_memory(ls); 
                MatrixFreeLinearSolver<Vector>::init_memory(ls); 
            }             
    };



    template<class Matrix, class Vector>
    class OperatorBasedQPSolver :   public virtual MatrixFreeQPSolver<Vector>,
                                    public virtual QPSolver<Matrix, Vector>
    {
    public:
        using MatrixFreeQPSolver<Vector>::update;
        using QPSolver<Matrix, Vector>::update;
        using MatrixFreeQPSolver<Vector>::solve;

        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;        

        virtual ~OperatorBasedQPSolver() {}

        virtual bool solve(const Matrix &A, const Vector &b, Vector &x) override
        {
            update(make_ref(A));
            return solve(operator_cast<Vector>(A), b, x);
        }

        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            QPSolver<Matrix, Vector>::update(op);
            update(operator_cast<Vector>(*op));
        }

        bool apply(const Vector &b, Vector &x) override
        {
            return solve(operator_cast<Vector>(*this->get_operator()), b, x);
        }

        virtual OperatorBasedQPSolver * clone() const override = 0;

        virtual void read(Input &in) override
        {
            MatrixFreeQPSolver<Vector>::read(in);
            QPSolver<Matrix, Vector>::read(in);
        }

        virtual void print_usage(std::ostream &os) const override
        {
            MatrixFreeQPSolver<Vector>::print_usage(os);
            QPSolver<Matrix, Vector>::print_usage(os);
        }

        virtual void init_memory(const SizeType & ls)  override 
        {
            MatrixFreeQPSolver<Vector>::init_memory(ls); 
            QPSolver<Matrix, Vector>::init_memory(ls); 
        }   

    };



}

#endif //UTOPIA_QP_SOLVER_HPP