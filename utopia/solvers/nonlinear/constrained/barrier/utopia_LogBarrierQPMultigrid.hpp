#ifndef UTOPIA_LOG_BARRIER_QP_MULTIGRID_SOLVER_HPP
#define UTOPIA_LOG_BARRIER_QP_MULTIGRID_SOLVER_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_PETSC

#include "utopia_BarrierMultigrid.hpp"
#include "utopia_Core.hpp"
#include "utopia_ILU.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_MatrixAgglomerator.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

#include "utopia_BlockAgglomerate.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class LogBarrierQPMultigrid final : public QPSolver<Matrix, Vector> {
        using Super = utopia::QPSolver<Matrix, Vector>;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using BarrierMultigrid = utopia::BarrierMultigrid<Matrix, Vector>;
        using QuadraticFunction = utopia::QuadraticFunction<Matrix, Vector>;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using MatrixAgglomerator = utopia::MatrixAgglomerator<Matrix>;

    public:
        LogBarrierQPMultigrid(const std::shared_ptr<LinearSolver> &linear_solver,
                              const std::shared_ptr<MatrixAgglomerator> &agglomerator = nullptr,
                              const std::shared_ptr<IterativeSolver> &linear_smoother = nullptr)
            : solver_(std::make_shared<BarrierMultigrid>(linear_solver)) {
            solver_->set_agglomerator(agglomerator);

            if (!linear_smoother) {
                solver_->set_linear_smoother(std::make_shared<ILU<Matrix, Vector>>());
            } else {
                solver_->set_linear_smoother(linear_smoother);
            }

            ensure_agglomerator();
            init_defaults();
        }

        LogBarrierQPMultigrid() {
            auto linear_solver = std::make_shared<OmniLinearSolver<Matrix, Vector>>();

            solver_ = std::make_shared<BarrierMultigrid>(linear_solver);

            auto linear_smoother = std::make_shared<ILU<Matrix, Vector>>();
            solver_->set_linear_smoother(linear_smoother);

            ensure_agglomerator();
            init_defaults();
        }

        LogBarrierQPMultigrid(const BarrierMultigrid &other) : Super(other) { solver_ = other.solver_->clone(); }

        ~LogBarrierQPMultigrid() {}
        LogBarrierQPMultigrid *clone() const override {
            auto ptr = utopia::make_unique<LogBarrierQPMultigrid>(*this);
            return ptr.release();
        }

        void read(Input &in) override {
            solver_->read(in);
            if (!Options()
                     .add_option("block_size", block_size_, "Used to specify the dimension of the vector problem")
                     .parse(in)) {
                return;
            }
        }

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

        void set_selection(const std::shared_ptr<Vector> &selection) override {
            assert(solver_);
            solver_->set_selection(selection);
        }

    private:
        std::shared_ptr<BarrierMultigrid> solver_;
        std::shared_ptr<QuadraticFunction> quadratic_function_;
        int block_size_{1};

        void init_defaults() {
            InputParameters params;
            params.set("use_coarse_space", true);
            params.set("debug", true);
            params.set("barrier_parameter", 1);
            params.set("barrier_parameter_shrinking_factor", 0.1);
            params.set("min_barrier_parameter", 1e-10);
            params.set("max_it", 40);
            params.set("barrier_function_type", "BoundedLogBarrier");

            params.set("use_non_linear_residual", false);
            params.set("pre_smoothing_steps", 5);
            params.set("post_smoothing_steps", 5);
            params.set("atol", 1e-7);
            params.set("rtol", 1e-8);
            params.set("mg_steps", 1);
            params.set("keep_initial_coarse_spaces", true);
            params.set("amg_n_coarse_spaces", 3);
            params.set("non_smooth_projection_weight", 0.9);

            solver_->read(params);
        }

        void ensure_agglomerator() {
            if (solver_->has_agglomerator()) return;

            switch (block_size_) {
                case 2: {
                    solver_->set_agglomerator(std::make_shared<BlockAgglomerate<Matrix, 2>>());
                    break;
                }
                case 3: {
                    solver_->set_agglomerator(std::make_shared<BlockAgglomerate<Matrix, 3>>());
                    break;
                }
                case 4: {
                    solver_->set_agglomerator(std::make_shared<BlockAgglomerate<Matrix, 4>>());
                    break;
                }
                default: {
                    solver_->set_agglomerator(std::make_shared<Agglomerate<Matrix>>());
                    break;
                }
            }
        }
    };  // namespace utopia

}  // namespace utopia
#endif  // UTOPIA_ENABLE_PETSC
#endif  // UTOPIA_LOG_BARRIER_QP_MULTIGRID_SOLVER_HPP
