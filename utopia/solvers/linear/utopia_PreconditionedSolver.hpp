#ifndef UTOPIA_PRECONDITIONED_SOLVER_HPP
#define UTOPIA_PRECONDITIONED_SOLVER_HPP

#include "utopia_IterativeSolver.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_PreconditionedSolverInterface.hpp"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class PreconditionedSolver : public IterativeSolver<Matrix, Vector>,
                                 virtual public PreconditionedSolverInterface<Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Preconditioner = utopia::Preconditioner<Vector>;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::PreconditionedSolverInterface<Vector> PreconditionedSolverInterface;

        using PreconditionedSolverInterface::update;

        virtual ~PreconditionedSolver() = default;

        // void update(const std::shared_ptr<Operator<Vector>> &op) override {
        // PreconditionedSolverInterface::update(op); }

        // For avoiding ambigous overrides with update(Operator)
        void update(const std::shared_ptr<Matrix> &op) { update(std::static_pointer_cast<const Matrix>(op)); }

        void update(const std::shared_ptr<const Operator<Vector>> &op) override {
            PreconditionedSolverInterface::update(op);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            IterativeSolver::update(op);
            if (this->precond_) {
                PreconditionedSolverInterface::update(std::static_pointer_cast<const Operator<Vector>>(op));
                auto ls_ptr = dynamic_cast<LinearSolver *>(this->precond_.get());
                if (ls_ptr) {
                    ls_ptr->update(op);
                }
            }
        }

        virtual void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec) {
            IterativeSolver::update(op);
            if (this->precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver *>(this->precond_.get());
                if (ls_ptr) {
                    ls_ptr->update(prec);
                } else {
                    auto delegate_ptr = dynamic_cast<DelegatePreconditioner<Matrix, Vector> *>(this->precond_.get());
                    if (delegate_ptr) {
                        delegate_ptr->update(prec);
                    }
                }
            }
        }

        void read(Input &in) override {
            IterativeSolver::read(in);
            PreconditionedSolverInterface::read(in);
        }

        void print_usage(std::ostream &os) const override {
            IterativeSolver::print_usage(os);
            PreconditionedSolverInterface::print_usage(os);
        }

        inline PreconditionedSolver &operator=(const PreconditionedSolver &other) {
            if (this == &other) {
                return *this;
            }

            IterativeSolver::operator=(other);
            PreconditionedSolverInterface::operator=(other);
            // This is/should called in the line above
            // this->copy_preconditioner_from(other);
            return *this;
        }

        PreconditionedSolver(const PreconditionedSolver &other)
            : PreconditionedSolverInterface(other), IterativeSolver(other) {
            // This is called in the constructor of PreconditionedSolverInterface
            // this->copy_preconditioner_from(other);
        }

        PreconditionedSolver() : PreconditionedSolverInterface(), IterativeSolver() {}

        PreconditionedSolver(PreconditionedSolver &&other) = delete;
    };
}  // namespace utopia

#endif  // UTOPIA_PRECONDITIONED_SOLVER_HPP
