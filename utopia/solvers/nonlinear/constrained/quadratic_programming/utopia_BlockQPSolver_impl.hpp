#ifndef UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP
#define UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP

#include <utility>

#include "utopia_Algorithms.hpp"
#include "utopia_BlockQPSolver.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::~BlockQPSolver() = default;

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::BlockQPSolver(std::shared_ptr<QPSolver> serial_solver)
        : serial_solver_(std::move(serial_solver)),
          local_lb_(std::make_shared<Vector>()),
          local_ub_(std::make_shared<Vector>()) {}

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::BlockQPSolver(const BlockQPSolver &other)
        : Super(other), local_lb_(std::make_shared<Vector>()), local_ub_(std::make_shared<Vector>()) {
        if (other.serial_solver_) {
            serial_solver_ = std::shared_ptr<QPSolver>(other.serial_solver_->clone());
        }
    }

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC> *BlockQPSolver<Matrix, Vector, PETSC>::clone() const {
        return new BlockQPSolver(*this);
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::read(Input &in) {
        QPSolver::read(in);

        if (serial_solver_) {
            // bool verb = serial_solver_->verbose();
            // serial_solver_->read(in);
            // serial_solver_->verbose(verb);
            in.get("block_solver", *serial_solver_);
        }
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::print_usage(std::ostream &os) const {
        QPSolver::print_usage(os);
    }

    template <class Matrix, class Vector>
    bool BlockQPSolver<Matrix, Vector, PETSC>::smooth(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BlockQPSolver::smooth");
        const Matrix &A = *this->get_operator();

        // TODO
        assert(false && "IMPLEMENT ME");

        return false;
    }

    template <class Matrix, class Vector>
    bool BlockQPSolver<Matrix, Vector, PETSC>::apply(const Vector &rhs, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BlockQPSolver::apply");

        if (rhs.comm().size() == 1) {
            serial_solver_->set_box_constraints(BoxConstraints<Vector>(this->lower_bound(), this->upper_bound()));
            return serial_solver_->apply(rhs, x);
        }

        this->init_solver("utopia BlockQPSolver comm.size = " + std::to_string(rhs.comm().size()),
                          {" it. ", "|| u - u_old ||"});
        // make sure we have both bounds
        this->fill_empty_bounds(layout(x));

        x_old_ = x;

        bool converged = false;

        int iteration = 0;
        while (!converged) {
            bool local_converged = step(rhs, x);
            if (!local_converged) {
                std::cerr << "[Warning] local stepper not converged" << std::endl;
            }

            c_ = x - x_old_;
            const Scalar diff = norm2(c_);

            if (this->verbose()) {
                PrintInfo::print_iter_status(iteration, {diff});
            }

            converged = this->check_convergence(iteration, 1, 1, diff);

            ++iteration;

            if (converged) break;

            x_old_ = x;
        }

        // x.comm().barrier();
        UTOPIA_TRACE_REGION_END("BlockQPSolver::apply");
        return converged;
    }

    template <class Matrix, class Vector>
    bool BlockQPSolver<Matrix, Vector, PETSC>::step(const Vector &rhs, Vector &x) {
        const Matrix &A = *this->get_operator();

        /////////////// global to local ///////////////

        // compute bound difference (ATTENTION !!!! residual_ is used as a temporary buffer)
        residual_ = this->get_lower_bound() - x;
        global_to_local(residual_, *local_lb_);

        residual_ = this->get_upper_bound() - x;
        global_to_local(this->get_upper_bound(), *local_ub_);

        // compute residual
        residual_ = A * x;
        residual_ = rhs - residual_;
        global_to_local(residual_, local_residual_);

        ////////////////// solve /////////////////////

        serial_solver_->set_box_constraints(BoxConstraints<Vector>(local_lb_, local_ub_));

        local_correction_.set(0.0);

        assert(local_residual_.comm().size() == 1);
        assert(local_correction_.comm().size() == 1);
        assert(local_lb_->comm().size() == 1);
        assert(local_ub_->comm().size() == 1);
        assert(serial_solver_->get_operator()->comm().size() == 1);

        bool converged = serial_solver_->apply(local_residual_, local_correction_);

        ///////////////// local to global /////////////////

        local_to_global_add(local_correction_, x);
        return converged;
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        local_correction_.zeros(serial_layout(layout.local_size()));
        assert(serial_solver_);
    }

    // template <class Matrix>
    // static void build_local_matrix(const Matrix &global, Matrix &local) {
    //     using SizeType = typename Traits<Matrix>::SizeType;
    //     using Scalar = typename Traits<Matrix>::Scalar;

    //     Matrix temp;
    //     local_block_view(global, temp);

    //     SizeType n_local = temp.rows();

    //     SizeType max_nnz = 0;
    //     std::vector<SizeType> nnz(n_local, 0);

    //     const SizeType offset = temp.row_range().begin();
    //     temp.read([&](const SizeType &i, const SizeType &, const Scalar &) { ++nnz[i - offset]; });

    //     if (!nnz.empty()) {
    //         max_nnz = *std::max_element(nnz.begin(), nnz.end());
    //     }

    //     local.sparse(serial_layout(n_local, n_local), max_nnz, 0);

    //     Write<Matrix> w_(local);
    //     temp.read([&](const SizeType &i, const SizeType &j, const Scalar &val) { local.set(i, j, val); });
    // }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op) {
        UTOPIA_TRACE_REGION_BEGIN("BlockQPSolver::update");
        IterativeSolver<Matrix, Vector>::update(op);

        if (op->comm().size() == 1) {
            serial_solver_->update(op);
            return;
        }

        std::shared_ptr<Matrix> A_view_ptr = std::make_shared<Matrix>();
        local_block_view(*op, *A_view_ptr);
        // build_local_matrix(*op, *A_view_ptr);

        serial_solver_->update(A_view_ptr);

        init_memory(row_layout(*op));
        UTOPIA_TRACE_REGION_END("BlockQPSolver::update");
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::global_to_local(const Vector &g, Vector &l) {
        if (empty(l)) {
            l.zeros(serial_layout(g.local_size()));
        }

        auto g_view = const_local_view_device(g);
        auto l_view = local_view_device(l);

        parallel_for(l.range_device(), UTOPIA_LAMBDA(const SizeType &i) { l_view.set(i, g_view.get(i)); });
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::local_to_global(const Vector &l, Vector &g) {
        if (empty(g)) {
            assert(false);
            m_utopia_error("g must be initialized outside");
            return;
        }

        auto g_view = local_view_device(g);
        auto l_view = const_local_view_device(l);

        parallel_for(l.range_device(), UTOPIA_LAMBDA(const SizeType &i) { g_view.set(i, l_view.get(i)); });
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::local_to_global_add(const Vector &l, Vector &g) {
        if (empty(g)) {
            assert(false);
            m_utopia_error("g must be initialized outside");
            return;
        }

        auto g_view = local_view_device(g);
        auto l_view = const_local_view_device(l);

        parallel_for(l.range_device(), UTOPIA_LAMBDA(const SizeType &i) { g_view.add(i, l_view.get(i)); });
    }

}  // namespace utopia

#endif  // UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP
