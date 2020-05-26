#ifndef UTOPIA_PETSC_EVAL_RESIDUAL_HPP
#define UTOPIA_PETSC_EVAL_RESIDUAL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template <class A, class X, class B>
    using PetscMatResidual = Binary<Tensor<B, 1>, Multiply<Tensor<A, 2>, Tensor<X, 1>>, Minus>;

    template <class A, class X, class B, class Traits>
    class Eval<PetscMatResidual<A, X, B>, Traits, PETSC> {
    public:
        typedef utopia::PetscMatResidual<A, X, B> Expr;
        using Result = X;

        inline static Result apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Result r;
            const auto &a = expr.right().left().derived();
            const auto &x = expr.right().right().derived();
            const auto &b = expr.left().derived();

            r.init(b.comm().get(), b.type(), b.local_size(), b.size());
            PetscErrorCode ierr;
            UTOPIA_UNUSED(ierr);
            ierr = MatResidual(a.raw_type(), b.raw_type(), x.raw_type(), r.raw_type());
            assert(ierr == 0);

            UTOPIA_TRACE_END(expr);
            return r;
        }
    };

    template <class Left, class A, class X, class B, class Traits>
    class Eval<Assign<Left, PetscMatResidual<A, X, B>>, Traits, PETSC> {
    public:
        typedef utopia::PetscMatResidual<A, X, B> Right;
        typedef utopia::Assign<Left, Right> Expr;

        inline static bool apply(const Expr &assign_expr) {
            UTOPIA_TRACE_BEGIN(assign_expr);
            auto &&expr = assign_expr.right();

            auto &r = assign_expr.left().derived();
            const auto &a = expr.right().left().derived();
            const auto &x = expr.right().right().derived();
            const auto &b = expr.left().derived();

            if (r.is_null() || r.size() != b.size()) {
                r.repurpose(b.comm().get(), b.type(), b.local_size(), b.size());
            }

            if (r.is_alias(b) || r.is_alias(x)) {
                auto temp = b;
                PetscErrorCode ierr;
                UTOPIA_UNUSED(ierr);
                ierr = MatResidual(a.raw_type(), b.raw_type(), x.raw_type(), temp.raw_type());
                assert(ierr == 0);
                r = std::move(temp);
            } else {
                PetscErrorCode ierr;
                UTOPIA_UNUSED(ierr);
                ierr = MatResidual(a.raw_type(), b.raw_type(), x.raw_type(), r.raw_type());
                assert(ierr == 0);
            }

            UTOPIA_TRACE_END(assign_expr);
            return true;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_EVAL_RESIDUAL_HPP
