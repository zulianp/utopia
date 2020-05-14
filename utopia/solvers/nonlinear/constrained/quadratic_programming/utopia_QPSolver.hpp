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
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        QPSolver() {}
        ~QPSolver() override {}
        QPSolver *clone() const override = 0;

        void init_memory(const Layout &layout) override {
            VariableBoundSolverInterface<Vector>::init_memory(layout);
            PreconditionedSolver<Matrix, Vector>::init_memory(layout);
        }
    };


    template<class Vector>
    class MatrixFreeQPSolver : public virtual MatrixFreeLinearSolver<Vector>,
                               public virtual VariableBoundSolverInterface<Vector> {
    public:
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        MatrixFreeQPSolver()
        {}

        ~MatrixFreeQPSolver() override {}

        MatrixFreeQPSolver *clone() const override = 0;

        void init_memory(const Layout &layout) override {
            VariableBoundSolverInterface<Vector>::init_memory(layout);
            MatrixFreeLinearSolver<Vector>::init_memory(layout);
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

        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        ~OperatorBasedQPSolver() override {}

        bool solve(const Matrix &A, const Vector &b, Vector &x) override {
            update(make_ref(A));
            return solve(operator_cast<Vector>(A), b, x);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            QPSolver<Matrix, Vector>::update(op);
            update(operator_cast<Vector>(*op));
        }

        bool apply(const Vector &b, Vector &x) override
        {
            return solve(operator_cast<Vector>(*this->get_operator()), b, x);
        }

        OperatorBasedQPSolver *clone() const override = 0;

        void read(Input &in) override {
            MatrixFreeQPSolver<Vector>::read(in);
            QPSolver<Matrix, Vector>::read(in);
        }

        void print_usage(std::ostream &os) const override {
            MatrixFreeQPSolver<Vector>::print_usage(os);
            QPSolver<Matrix, Vector>::print_usage(os);
        }

        void init_memory(const Layout &layout) override {
            MatrixFreeQPSolver<Vector>::init_memory(layout);
            QPSolver<Matrix, Vector>::init_memory(layout);
        }
    };

}

#endif //UTOPIA_QP_SOLVER_HPP