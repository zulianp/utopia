#ifndef UTOPIA_PETSC_BDD_QP_SOLVER_HPP
#define UTOPIA_PETSC_BDD_QP_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_QPSolver.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class BDDQPSolver final : public QPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using MatrixFreeQPSolver = utopia::MatrixFreeQPSolver<Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

    public:
        void set_solver(const std::shared_ptr<MatrixFreeQPSolver> &solver);

        explicit BDDQPSolver(std::shared_ptr<MatrixFreeQPSolver> solver);

        BDDQPSolver();

        ~BDDQPSolver() override;
        BDDQPSolver *clone() const override;

        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;

        bool smooth(const Vector &, Vector &) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;

        void set_selection(const std::shared_ptr<Vector> &boolean_selector) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        void determine_boolean_selector();
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_BDD_QP_SOLVER_HPP
