#ifndef UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_HPP
#define UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_QPSolver.hpp"

#include "utopia_polymorphic_LinearSolver.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class PrimalInteriorPointSolver final : public QPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

    public:
        inline void set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver) {
            linear_solver_ = linear_solver;
        }

        explicit PrimalInteriorPointSolver(
            std::shared_ptr<LinearSolver> linear_solver = std::make_shared<OmniLinearSolver<Matrix, Vector>>());
        ~PrimalInteriorPointSolver() override;
        PrimalInteriorPointSolver *clone() const override;

        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;

        bool smooth(const Vector &, Vector &) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        void fallback_solver(const std::shared_ptr<LinearSolver> &fallback_solver);

    private:
        class Impl;

        std::shared_ptr<LinearSolver> linear_solver_;
        std::unique_ptr<Impl> details_;
    };

}  // namespace utopia

#endif  // UTOPIA_PRIMAL_INTERIOR_POINT_SOLVER_HPP