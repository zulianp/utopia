#ifndef UTOPIA_BLOCK_QP_SOLVER_HPP
#define UTOPIA_BLOCK_QP_SOLVER_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_SolverForwardDeclarations.hpp"

#include <cmath>
#include "utopia_BlockQPSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class BlockQPSolver {};

    // FIXME make it work for all other backends
    template <class Matrix, class Vector>
    class BlockQPSolver<Matrix, Vector, PETSC> : public QPSolver<Matrix, Vector> {
       public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using Comm = typename Traits<Vector>::Communicator;
        using QPSolver = utopia::QPSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

        ~BlockQPSolver();
        BlockQPSolver(const std::shared_ptr<QPSolver> &serial_solver);
        BlockQPSolver(const BlockQPSolver &other);

        virtual BlockQPSolver *clone() const override;
        virtual void read(Input &in) override;
        virtual void print_usage(std::ostream &os) const override;
        virtual bool smooth(const Vector &b, Vector &x) override;
        virtual bool apply(const Vector &b, Vector &x) override;
        virtual void init_memory(const Layout &layout) override;
        virtual void update(const std::shared_ptr<const Matrix> &op) override;

       private:
        std::shared_ptr<QPSolver> serial_solver_;
        Vector local_rhs_, local_x_;

        static void global_to_local(const Vector &g, Vector &l);
        static void local_to_global(const Vector &l, Vector &g);
    };
}  // namespace utopia

#endif  // UTOPIA_BLOCK_QP_SOLVER_HPP
