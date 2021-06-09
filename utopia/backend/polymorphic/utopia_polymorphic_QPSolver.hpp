#ifndef UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_QPSolver.hpp"

#include <memory>

namespace utopia {

    template <class Matrix, class Vector>
    class OmniQPSolver : public QPSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::QPSolver<Matrix, Vector> Super;
        using BoxConstraints = utopia::BoxConstraints<Vector>;

    public:
        OmniQPSolver();
        ~OmniQPSolver() override;
        OmniQPSolver *clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void read(Input &in) override;

    private:
        std::unique_ptr<QPSolver<Matrix, Vector>> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
