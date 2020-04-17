#ifndef UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP
#define UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_BlockQPSolver.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::~BlockQPSolver() {}

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::BlockQPSolver(const std::shared_ptr<QPSolver> &serial_solver)
        : serial_solver_(serial_solver) {}

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC>::BlockQPSolver(const BlockQPSolver &other) {
        if (other.serial_solver_) {
            serial_solver_ = std::shared_ptr<QPSolver>(other.serial_solver_->clone());
        }
    }

    template <class Matrix, class Vector>
    BlockQPSolver<Matrix, Vector, PETSC> *BlockQPSolver<Matrix, Vector, PETSC>::clone() const {
        auto ptr = new BlockQPSolver(*this);
        ptr->set_box_constraints(this->get_box_constraints());
        return ptr;
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::read(Input &in) {
        QPSolver::read(in);

        if (serial_solver_) {
            serial_solver_->read(in);
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
    bool BlockQPSolver<Matrix, Vector, PETSC>::apply(const Vector &b, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BlockQPSolver::apply");

        this->init_solver("utopia BlockQPSolver", {" it. ", "|| u - u_old ||"});

        const Matrix &A = *this->get_operator();

        bool converged = false;

        // global to local

        // solve

        // local to global

        // TODO
        UTOPIA_TRACE_REGION_END("BlockQPSolver::apply");
        return converged;
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::init_memory(const Layout &layout) {
        Super::init_memory(layout);
        assert(serial_solver_);
    }

    template <class Matrix, class Vector>
    void BlockQPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op) {
        UTOPIA_TRACE_REGION_BEGIN("BlockQPSolver::update");
        IterativeSolver<Matrix, Vector>::update(op);

        // TODO
        std::shared_ptr<Matrix> A_view_ptr = std::make_shared<Matrix>();
        local_block_view(*op, *A_view_ptr);

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

}  // namespace utopia

#endif  // UTOPIA_BLOCK_QP_SOLVER_IMPL_HPP
