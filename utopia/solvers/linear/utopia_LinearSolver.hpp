#ifndef UTOPIA_SOLVER_SOLVER_H
#define UTOPIA_SOLVER_SOLVER_H

#include <cassert>
#include <string>

#include "utopia_Clonable.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Core.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Utils.hpp"

namespace utopia {
    /**
     * @brief      The base class for linear solvers.
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class LinearSolver : virtual public Preconditioner<Vector> {
    public:
        using Super = utopia::Preconditioner<Vector>;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        using Super::update;

        LinearSolver() = default;
        ~LinearSolver() override = default;
        LinearSolver &operator=(const LinearSolver &) = default;
        LinearSolver &operator=(LinearSolver &&) = default;
        LinearSolver(LinearSolver &&other) = default;

        LinearSolver(const LinearSolver &other) {
            UTOPIA_UNUSED(other);
            assert((!other.op_) && "cannot be copied once initialized");
        }

        bool apply(const Vector &rhs, Vector &sol) override = 0;

        void read(Input &in) override { Preconditioner<Vector>::read(in); }

        void print_usage(std::ostream &os) const override { Preconditioner<Vector>::print_usage(os); }

        /**
         * @brief      Solve routine.
         * @param[in]  A
         * @param[in]  b
         * @param      x0
         *
         * @return
         */
        virtual bool solve(const Matrix &A, const Vector &b, Vector &x0) {
            update(make_ref(A));
            return apply(b, x0);
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        virtual void update(const std::shared_ptr<const Matrix> &op)  // override
        {
            op_ = op;
        }

        inline const std::shared_ptr<const Matrix> &get_operator() const {
            assert(op_);
            return op_;
        }

        inline bool has_operator() const { return static_cast<bool>(op_); }

        LinearSolver *clone() const override = 0;

    private:
        std::shared_ptr<const Matrix> op_;
    };

}  // namespace utopia

#endif  // UTOPIA_SOLVER_SOLVER_H
