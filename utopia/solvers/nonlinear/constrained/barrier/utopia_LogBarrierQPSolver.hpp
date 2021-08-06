#ifndef UTOPIA_LOG_BARRIER_QP_SOLVER_HPP
#define UTOPIA_LOG_BARRIER_QP_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_LogBarrierSolver.hpp"
#include "utopia_Newton.hpp"
#include "utopia_QuadraticFunction.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class LogBarrierQPSolver final : public QPSolver<Matrix, Vector> {
        using Super = utopia::QPSolver<Matrix, Vector>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LogBarrierSolver = utopia::LogBarrierSolver<Matrix, Vector>;
        using QuadraticFunction = utopia::QuadraticFunction<Matrix, Vector>;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;

    public:
        LogBarrierQPSolver(const std::shared_ptr<LinearSolver> &linear_solver =
                               std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>())
            : solver_(std::make_shared<LogBarrierSolver>(linear_solver)) {}

        LogBarrierQPSolver(const LogBarrierSolver &other) : Super(other) { solver_ = other.solver_->clone(); }

        ~LogBarrierQPSolver() {}
        LogBarrierQPSolver *clone() const override {
            auto ptr = utopia::make_unique<LogBarrierQPSolver>(*this);
            return ptr.release();
        }

        void read(Input &in) override { solver_->read(in); }
        void print_usage(std::ostream &os) const override { solver_->print_usage(os); }

        bool smooth(const Vector &, Vector &) override {
            assert(false && "IMPLEMENT ME");
            return false;
        }

        bool apply(const Vector &b, Vector &x) override {
            quadratic_function_ = std::make_shared<QuadraticFunction>(this->get_operator(), make_ref(b));
            solver_->set_box_constraints(this->get_box_constraints());
            return solver_->solve(*quadratic_function_, x);
        }

        void init_memory(const Layout &layout) override { Super::init_memory(layout); }
        void update(const std::shared_ptr<const Matrix> &op) override { Super::update(op); }

    private:
        std::shared_ptr<LogBarrierSolver> solver_;
        std::shared_ptr<QuadraticFunction> quadratic_function_;
    };

}  // namespace utopia

#endif  // UTOPIA_LOG_BARRIER_QP_SOLVER_HPP
