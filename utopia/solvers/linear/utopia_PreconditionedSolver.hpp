#ifndef UTOPIA_PRECONDITIONED_SOLVER_HPP
#define UTOPIA_PRECONDITIONED_SOLVER_HPP

#include "utopia_IterativeSolver.hpp"
#include "utopia_LinearSolver.hpp"

namespace utopia {
    template<class Matrix, class Vector>
    class PreconditionedSolver : public IterativeSolver<Matrix, Vector> {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

        /**
        * @brief      Sets the preconditioner.
        *
        * @param[in]  precond  The precondition
        */
        virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
        {
            precond_ = precond;
        }

         /**
         * @brief      Resets the preconditioner.
         */
        void reset_preconditioner()
        {
            precond_.reset();
        }

        const std::shared_ptr<Preconditioner> &get_preconditioner() const
        {
            return precond_;
        }

        inline bool has_preconditioner() const
        {
            return static_cast<bool>(precond_);
        }

        void update(const std::shared_ptr<const Matrix> &op) override { update(op, op); }

        virtual void update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec)
        {
            IterativeSolver::update(op);

            if(precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver *>(precond_.get());
                if(ls_ptr) {
                    ls_ptr->update(prec);
                } else {
                    auto delegate_ptr = dynamic_cast< DelegatePreconditioner<Matrix, Vector> *>(precond_.get());
                    if(delegate_ptr) {
                        delegate_ptr->update(prec);
                    }
                }
            }
        }

        void read(Input &in) override {
            IterativeSolver::read(in);
            if(precond_) {
                in.get("precond", *precond_);
            }
        }

        void print_usage(std::ostream &os) const override {
            IterativeSolver::print_usage(os);
            this->print_param_usage(os, "precond", "Preconditioner", "Input parameters for preconditioner.", "-");
        }

        inline PreconditionedSolver &operator=(const PreconditionedSolver &other)
        {
            if(this == &other) {
                return *this;
            }

            IterativeSolver::operator=(other);
            copy_preconditioner_from(other);
            return *this;
        }

        PreconditionedSolver(const PreconditionedSolver &other)
        : IterativeSolver(other)
        {
            copy_preconditioner_from(other);
        }

        PreconditionedSolver() : IterativeSolver() {}

    private:
        std::shared_ptr<Preconditioner> precond_;

        void copy_preconditioner_from(const PreconditionedSolver &other)
        {
            if(other.has_preconditioner()) {
                this->set_preconditioner(std::shared_ptr<Preconditioner>(other.precond_->clone()));
            }
        }

    };
}

#endif //UTOPIA_PRECONDITIONED_SOLVER_HPP
