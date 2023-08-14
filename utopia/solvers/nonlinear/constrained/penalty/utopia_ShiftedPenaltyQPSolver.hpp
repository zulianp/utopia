#ifndef UTOPIA_SHIFTED_PENALTY_QP_SOLVER_HPP
#define UTOPIA_SHIFTED_PENALTY_QP_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_QPSolver.hpp"

#include <memory>
#include <vector>

namespace utopia {

    // Method from DOI 10.1007/s00466-015-1150-5
    template <class Matrix, class Vector = typename Traits<Matrix>::Vector, int Backend = Traits<Vector>::Backend>
    class ShiftedPenaltyQPSolver final : public QPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

    public:
        void set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver);

        explicit ShiftedPenaltyQPSolver(std::shared_ptr<LinearSolver> linear_solver);

        ShiftedPenaltyQPSolver();

        ~ShiftedPenaltyQPSolver() override;
        ShiftedPenaltyQPSolver *clone() const override;

        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;

        bool smooth(const Vector &, Vector &) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;

        void set_selection(const std::shared_ptr<Vector> &boolean_selector) override;

        UTOPIA_NVCC_PRIVATE
        class Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_SHIFTED_PENALTY_QP_SOLVER_HPP
