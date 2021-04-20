#ifndef UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP

#include <memory>
#include "utopia_Describable.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class OmniLinearSolver : public LinearSolver<Matrix, Vector>, public Describable {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Super = utopia::LinearSolver<Matrix, Vector>;
        using Super::update;

        inline void describe(std::ostream &os) const override { os << "OmniLinearSolver"; }

    public:
        OmniLinearSolver();
        ~OmniLinearSolver() override;
        OmniLinearSolver *clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void update(const std::shared_ptr<const Matrix> &mat) override;
        void read(Input &in) override;

    private:
        std::unique_ptr<LinearSolver<Matrix, Vector>> impl_;
    };

}  // namespace utopia
#endif  // UTOPIA_POLYMORPHIC_LINEAR_SOLVER_HPP
