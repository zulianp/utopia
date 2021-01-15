#ifndef UTOPIA_UTOPIA_PRECONDITIONER_HPP
#define UTOPIA_UTOPIA_PRECONDITIONER_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Expression.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Input.hpp"
#include "utopia_Memory.hpp"
#include "utopia_Operator.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_make_unique.hpp"

#include <cassert>
#include <memory>
#include <utility>

// #define UTOPIA_W_VECTOR(Tensor) utopia::Wrapper<typename
// utopia::Traits<Tensor>::Vector, 1>

namespace utopia {

    template <class Vector, class Fun>
    class LambdaOperator final : public Operator<Vector> {
    public:
        using Communicator = typename Traits<Vector>::Communicator;

        LambdaOperator(Communicator &comm, Size size, Size local_size, Fun fun)
            : comm_(comm), size_(std::move(size)), local_size_(std::move(local_size)), fun_(fun) {}

        inline bool apply(const Vector &rhs, Vector &sol) const override { return fun_(rhs, sol); }

        inline Size size() const override { return size_; }

        inline Size local_size() const override { return local_size_; }

        inline Communicator &comm() override { return comm_; }

        inline const Communicator &comm() const override { return comm_; }

    private:
        Communicator &comm_;
        Size size_, local_size_;
        Fun fun_;
    };

    template <typename Vector>
    std::unique_ptr<LambdaOperator<Vector, std::function<bool(const Vector &, Vector &)>>> op(
        typename Traits<Vector>::Communicator &comm,
        const Size &size,
        const Size &local_size,
        std::function<bool(const Vector &, Vector &)> f) {
        return utopia::make_unique<LambdaOperator<Vector, std::function<bool(const Vector &, Vector &)>>>(
            comm, size, local_size, f);
    }

    template <class Vector>
    class Preconditioner : public virtual Configurable,
                           public virtual Clonable,
                           public virtual MemoryInterface<Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        ~Preconditioner() override = default;
        virtual bool apply(const Vector &rhs, Vector &sol) = 0;

        void read(Input & /*in*/) override {}

        void print_usage(std::ostream & /*os*/) const override {}

        // this should be left as virtual once all solvers allocations improve
        void init_memory(const Layout & /*ls*/) override {}

        // TODO
        virtual void update(const Operator<Vector> &A) { UTOPIA_UNUSED(A); }

        Preconditioner *clone() const override = 0;
    };

    template <class Expr, class Vector>
    class ExprPreconditioner : public Preconditioner<Vector> {
    public:
        bool apply(const Vector &rhs, Vector &sol) override {
            sol = expr_ * rhs;
            return true;
        }

        ExprPreconditioner(const Expr &expr) : expr_(expr) {}

        ExprPreconditioner *clone() const override { return new ExprPreconditioner(*this); }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template <class Vector, class Derived>
    std::shared_ptr<ExprPreconditioner<Derived, Vector>> make_preconditioner(const Expression<Derived> &expr) {
        return std::make_shared<ExprPreconditioner<Derived, Vector>>(expr.derived());
    }

    template <class Matrix, class Vector>
    class DelegatePreconditioner : public Preconditioner<Vector> {
    public:
        using Preconditioner<Vector>::update;

        bool apply(const Vector & /*rhs*/, Vector & /*sol*/) override {
            std::cerr << "[Warning] DelegatePreconditioner::apply doing nothing ... \n";
            return true;
        }

        void update(const std::shared_ptr<const Matrix> &op) { op_ = op; }

        const std::shared_ptr<const Matrix> &get_matrix() const { return op_; }

        DelegatePreconditioner *clone() const override { return new DelegatePreconditioner(*this); }

    private:
        std::shared_ptr<const Matrix> op_;
    };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_PRECONDITIONER_HPP
