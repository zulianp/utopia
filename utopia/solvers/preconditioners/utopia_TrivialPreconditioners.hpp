#ifndef UTOPIA_TRIVIAL_PRECONDITIONERS_HPP
#define UTOPIA_TRIVIAL_PRECONDITIONERS_HPP

#include "utopia_Allocations.hpp"
#include "utopia_SolverForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

#include <cassert>
#include <memory>

#define UTOPIA_W_VECTOR(Tensor) utopia::Wrapper<typename utopia::Traits<Tensor>::Vector, 1>

namespace utopia {

    template <class Matrix, class Vector>
    class InvDiagPreconditioner final : public LinearSolver<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        bool apply(const Vector &rhs, Vector &sol) override {
            UTOPIA_NO_ALLOC_BEGIN("InvDiagPreconditioner:region3");
            sol = e_mul(d, rhs);
            UTOPIA_NO_ALLOC_END();
            return true;
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override {
            LinearSolver<Matrix, Vector>::update(op);
            UTOPIA_NO_ALLOC_BEGIN("InvDiagPreconditioner:region1");
            d = diag(*op);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("InvDiagPreconditioner:region2");
            // d = 1.0 / d;

            e_pseudo_inv(d, d);
            UTOPIA_NO_ALLOC_END();
        }

        Vector get_d() { return d; }

        InvDiagPreconditioner *clone() const override { return new InvDiagPreconditioner(*this); }

        void init_memory(const Layout &layout) override { d.zeros(layout); }

    private:
        Vector d;
    };

    template <class Vector>
    class IdentityPreconditioner final : public Preconditioner<Vector> {
    public:
        bool apply(const Vector &rhs, Vector &sol) override {
            sol = rhs;
            return true;
        }

        IdentityPreconditioner *clone() const override { return new IdentityPreconditioner(*this); }
    };

    template <class Vector>
    class VectorPreconditioner final : public Preconditioner<Vector> {
    public:
        bool apply(const Vector &rhs, Vector &sol) override {
            sol = e_mul(*vec_, rhs);
            return true;
        }

        VectorPreconditioner *clone() const override { return new VectorPreconditioner(*this); }

        VectorPreconditioner(const std::shared_ptr<Vector> &vec) : vec_(vec) {}

    private:
        std::shared_ptr<Vector> vec_;
    };

    template <class Vector>
    class FunctionPreconditioner final : public Preconditioner<Vector> {
    public:
        FunctionPreconditioner(const std::function<void(const Vector &, Vector &)> operator_action)
            : operator_action_(operator_action) {}

        bool apply(const Vector &rhs, Vector &ret) override {
            operator_action_(rhs, ret);
            return true;
        }

        FunctionPreconditioner *clone() const override { return new FunctionPreconditioner(*this); }

    private:
        std::function<void(const Vector &, Vector &)> operator_action_;
    };

}  // namespace utopia

#endif  // UTOPIA_TRIVIAL_PRECONDITIONERS_HPP