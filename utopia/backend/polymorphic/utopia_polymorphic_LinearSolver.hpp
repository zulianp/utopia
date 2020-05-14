#ifndef UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP

#include "utopia_LinearSolver.hpp"
#include "utopia_Traits.hpp"
#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class PolymorphicLinearSolver : public LinearSolver<Matrix, Vector>  {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Super;

    public:

        PolymorphicLinearSolver();
        ~PolymorphicLinearSolver() override;
        PolymorphicLinearSolver * clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void update(const std::shared_ptr<const Matrix> &mat) override;
        void read(Input &in) override;

    private:
        std::unique_ptr<LinearSolver<Matrix, Vector>> impl_;

    };

}
#endif //UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP
