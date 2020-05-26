#ifndef UTOPIA_ABSTRACT_LINEAR_SOLVER_HPP
#define UTOPIA_ABSTRACT_LINEAR_SOLVER_HPP

#include "utopia_AbstractMatrix.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_LinearSolver.hpp"

namespace utopia {

    template <typename Scalar, typename SizeType>
    class AbstractLinearSolver
        : public LinearSolver<AbstractMatrix<Scalar, SizeType>, AbstractVector<Scalar, SizeType>> {
    public:
        ~AbstractLinearSolver() override = default;
    };

    template <class Solver>
    class LinearSolverWrapper
        : public AbstractLinearSolver<typename Traits<Solver>::Scalar, typename Traits<Solver>::SizeType> {
    public:
        using Matrix = typename Traits<Solver>::Matrix;
        using Vector = typename Traits<Solver>::Vector;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;
        using AbstractMatrix = utopia::AbstractMatrix<Scalar, SizeType>;
        using Super = utopia::AbstractLinearSolver<Scalar, SizeType>;

        LinearSolverWrapper() : impl_(std::make_shared<Solver>()) {}

        LinearSolverWrapper(Solver *solver_ptr) : impl_(std::shared_ptr<Solver>(solver_ptr)) {}

        // template<class... Args>
        // LinearSolverWrapper(Args &&...args)
        // : impl_(std::make_shared<Solver>(std::forward<Args>(args)...))
        // {}

        template <class... Args>
        void construct(Args &&... args) {
            impl_ = std::make_shared<Solver>(std::forward<Args>(args)...);
        }

        inline bool apply(const AbstractVector &rhs, AbstractVector &sol) override {
            return impl_->apply(extract_vector(rhs), extract_vector(sol));
        }

        inline void read(Input &in) override { impl_->read(in); }

        inline void print_usage(std::ostream &os) const override { impl_->print_usage(os); }

        /**
         * @brief      Solve routine.
         * @param[in]  A
         * @param[in]  b
         * @param      x0
         *
         * @return
         */
        inline bool solve(const AbstractMatrix &A, const AbstractVector &rhs, AbstractVector &sol) override {
            update(make_ref(A));
            return impl_->apply(extract_vector(rhs), extract_vector(sol));
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        inline void update(const std::shared_ptr<const AbstractMatrix> &op) override {
            auto w = std::dynamic_pointer_cast<const Wrapper<Matrix>>(op);
            Super::update(w);
            impl_->update(w->ptr());
        }

        inline LinearSolverWrapper *clone() const override {
            auto ptr = utopia::make_unique<LinearSolverWrapper>(impl_->clone());
            return ptr.release();
        }

    private:
        std::shared_ptr<Solver> impl_;

        inline Vector &extract_vector(AbstractVector &v) const {
            auto w = dynamic_cast<Wrapper<Vector> *>(&v);
            // FIXME (use eventually local buffer)
            assert(w);
            return w->get();
        }

        inline const Vector &extract_vector(const AbstractVector &v) const {
            auto w = dynamic_cast<const Wrapper<Vector> *>(&v);
            // FIXME (use eventually local buffer)
            assert(w);
            return w->get();
        }
    };

}  // namespace utopia

#endif  // UTOPIA_ABSTRACT_LINEAR_SOLVER_HPP
