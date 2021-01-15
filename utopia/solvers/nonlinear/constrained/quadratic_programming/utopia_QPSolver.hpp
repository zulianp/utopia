#ifndef UTOPIA_QP_SOLVER_HPP
#define UTOPIA_QP_SOLVER_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"

#include <cmath>

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolver : public virtual VariableBoundSolverInterface<Vector>, public PreconditionedSolver<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        QPSolver(const QPSolver &other)
            : VariableBoundSolverInterface<Vector>(other), PreconditionedSolver<Matrix, Vector>(other) {}

        QPSolver() = default;
        ~QPSolver() override = default;
        QPSolver *clone() const override = 0;

        void init_memory(const Layout &layout) override {
            VariableBoundSolverInterface<Vector>::init_memory(layout);
            PreconditionedSolver<Matrix, Vector>::init_memory(layout);
        }
    };

    template <class Vector>
    class MatrixFreeQPSolver : public virtual VariableBoundSolverInterface<Vector>,
                               public MatrixFreeLinearSolver<Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        MatrixFreeQPSolver(const MatrixFreeQPSolver &other)
            : VariableBoundSolverInterface<Vector>(other), MatrixFreeLinearSolver<Vector>(other) {}

        MatrixFreeQPSolver() = default;

        ~MatrixFreeQPSolver() override = default;

        MatrixFreeQPSolver *clone() const override = 0;

        void init_memory(const Layout &layout) override {
            VariableBoundSolverInterface<Vector>::init_memory(layout);
            MatrixFreeLinearSolver<Vector>::init_memory(layout);
        }
    };

    template <class Matrix, class Vector>
    class OperatorBasedQPSolver : public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector> {
    public:
        using MatrixFreeQPSolver<Vector>::update;
        using QPSolver<Matrix, Vector>::update;
        using MatrixFreeQPSolver<Vector>::solve;

        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        OperatorBasedQPSolver() = default;

        // FIX since VariableBoundSolverInterface<Vector>(other) creates problems for
        // NVCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreorder"
        OperatorBasedQPSolver(const OperatorBasedQPSolver &other)
            : VariableBoundSolverInterface<Vector>(other),
              MatrixFreeQPSolver<Vector>(other),
              QPSolver<Matrix, Vector>(other) {}

#pragma GCC diagnostic pop

        OperatorBasedQPSolver &operator=(const OperatorBasedQPSolver &other) {
            if (this == &other) return *this;

            VariableBoundSolverInterface<Vector>::operator=(other);
            MatrixFreeQPSolver<Vector>::operator=(other);
            QPSolver<Matrix, Vector>::operator=(other);

            return *this;
        }

        ~OperatorBasedQPSolver() override = default;

        bool solve(const Matrix &A, const Vector &b, Vector &x) override {
            update(make_ref(A));
            return solve(operator_cast<Vector>(A), b, x);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            QPSolver<Matrix, Vector>::update(op);
            update(operator_cast<Vector>(*op));
        }

        bool apply(const Vector &b, Vector &x) override {
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

}  // namespace utopia

#endif  // UTOPIA_QP_SOLVER_HPP
